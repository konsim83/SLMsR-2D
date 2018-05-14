function solve_FEM_periodic_square(par :: Parameter.Parameter_FEM, problem :: T) where {T<:Problem.AbstractPhysicalProblem}

    # Build mesh of unit square (0,1)x(0,1)
    mesh = Mesh.mesh_unit_square(par.n_edge_per_seg)
    if par.n_refinement>0
        mesh = Mesh.refine_rg(mesh, par.n_refinement)

        # Subdividing edges messes up boundary markers. We need to correct
        # that.
        ind_point_boundary = sort(unique(mesh.edge[:,mesh.edge_marker.!=0]))
        mesh.point_marker[:] = zeros(Int, size(mesh.point_marker))
        mesh.point_marker[ind_point_boundary] = ones(Int, size(ind_point_boundary))
    end

    # Set up reference element
    ref_el = FEM.RefEl_Pk(par.n_order_FEM)

    # Set up quadrature rule
    quad = Quad.Quad_simplex(par.n_order_quad)

    # Set up degrees of freedom handler
    periodicityInfo = Mesh.DoublePeriodicUnitSquare()
    dof = FEM.Dof_Pk_periodic(mesh, problem, periodicityInfo, ref_el.n_order)

    # Set up solution structure
    solution = FEM.Solution_FEM(dof, par)

    # Set up time integrator
    timeStepper = TimeIntegrator.ImplEuler(dof, mesh, problem)
    # +++++++++++++++++++

    # Call actual solver. Pass solution data by reference.       
    solve_problem_FEM!(mesh,
                       ref_el,
                       dof,
                       quad,
                       timeStepper,
                       par,
                       problem,
                       solution)
    
    return solution, mesh
end


# ------------------------------------------------------------------------------------------------------------------------------


function solve_problem_FEM!(mesh :: Mesh.TriangleMesh.TriMesh,
                            ref_el :: FEM.AbstractRefEl,
                            dof :: FEM.AbstractDof,
                            quad :: Quad.AbstractQuad,
                            timeStepper :: TimeIntegrator.AbstractTimeIntegrator,
                            par :: Parameter.AbstractParameter,
                            problem :: Problem.AbstractProblem,
                            solution :: FEM.Solution_FEM)

    println("\n\n--------------------------------------------------------------")
    println("Computing standard FEM solution:\n")
    
    println("\t Mesh type:   $(mesh.mesh_info)")
    if dof.is_periodic
        println("\t Problem type:   periodic --- $(problem.type_info)")
    else
        println("\t Problem type:   non-periodic --- $(problem.type_info)")
    end
    println("\t Time intergrator:   $(typeof(timeStepper))")
    println("\t number of elements:   $(mesh.n_cell)")
    println("\t number of time steps:   $(par.n_steps)")
    println("\t Number of active dofs:   $(dof.n_true_dof - dof.n_node_dirichlet)")
    println("\t Number of constraint dofs:   $(dof.n_node_dirichlet)")
    
    # Setup initial data
    solution.u[:,1] = Problem.u_init(problem, mesh.point)

    # Make step from k_time to k_time+1
    M, A, D = [], [] ,[]
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    for k_time in 1:par.n_steps
        # Current time index
        Time = k_time*par.dt


        # Assemble the system matrices only when necessary
        if k_time==1
            M = FEM.assemble_mass(mesh,
                                   dof,
                                   ref_el,
                                   quad,
                                   problem)
        end
        
        if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
            A = FEM.assemble_advection(mesh,
                                         dof,
                                         ref_el,
                                         quad,
                                         problem,
                                         Time+par.dt)
        end

        if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
            D = FEM.assemble_diffusion(mesh,
                                         dof,
                                         ref_el,
                                         quad,
                                         problem,
                                         Time+par.dt)
        end

        # Zero forcing
        f = zeros(mesh.n_point)

        # Set the system matrices 
        TimeIntegrator.updateSystem!(timeStepper.systemData, M, D-A, f, dof, solution.u[:,k_time], par.dt)

        # Make a single time step
        TimeIntegrator.makeStep!(timeStepper, dof, view(solution.u, :, k_time+1), solution.u[:,k_time])

        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")

    return nothing
end