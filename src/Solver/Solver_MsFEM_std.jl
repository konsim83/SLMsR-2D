function solve_MsFEM_periodic_square(par :: Parameter.Parameter_MsFEM, problem :: T) where {T<:Problem.AbstractProblem}
    
    # Build mesh of unit square (0,1)x(0,1)
    m_coarse = Mesh.mesh_unit_square(par.n_edge_per_seg)
    m_simplex = Mesh.mesh_unit_simplex()
    if par.n_refinement>0
        # If this is not the case then the conformal MsFEM is just the
        # standard FEM
        m_simplex = Mesh.refine_rg(m_simplex, par.n_refinement)

        # Subdividing edges messes up boundary markers. We need to correct
        # that.
        ind_point_boundary = sort(unique(m_simplex.edge[:,m_simplex.edge_marker.!=0]))
        m_simplex.point_marker[:] = zeros(Int, size(m_simplex.point_marker))
        m_simplex.point_marker[ind_point_boundary] = ones(Int, size(ind_point_boundary))
    end
    
    mesh_collection = Mesh.TriMesh_collection(m_coarse, m_simplex)
    #mesh_collection = Mesh.TriMesh_collection(m_coarse, par.n_edge_per_seg_f)
    

    # Set up reference element
    ref_el = FEM.RefEl_Pk(1)
    ref_el_f = FEM.RefEl_Pk(par.n_order_FEM_f)

    # Set up quadrature rule
    quad_f = Quad.Quad_simplex(par.n_order_quad_f)

    # Set up degrees of freedom handler
    periodicityInfo = Mesh.DoublePeriodicUnitSquare()
    dof = FEM.Dof_Pk_periodic(mesh_collection.mesh, problem, periodicityInfo, 1)
    problem_f = Array{Problem.AbstractBasisProblem,1}(undef, mesh_collection.mesh.n_cell)
    for i_cell in 1:mesh_collection.mesh.n_cell
        # Set up local problem by geometry
        tri = Geometry.Triangle(FEM.map_ref_point(dof, ref_el.node, i_cell))
        problem_f[i_cell] = Problem.BasisFun(problem, tri)
    end
    dof_collection = FEM.Dof_collection(mesh_collection, dof, problem_f, ref_el_f.n_order)

    # Set up solution structure
    solution = FEM.Solution_MsFEM(dof_collection, par)


    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # Solve the local problems
    for i_cell in 1:mesh_collection.mesh.n_cell        
        # Set up time integrator
        timeStepper = TimeIntegrator.ImplEuler(dof_collection.dof_f[i_cell],
                                                mesh_collection.mesh_f[i_cell],
                                                problem_f[i_cell])
        # Call actual solver. Pass solution data by reference.  Solves
        # the i-th cell problem.
        solve_problem_MsFEM_basis!(mesh_collection.mesh_f[i_cell],
                                     ref_el_f,
                                     dof_collection.dof_f[i_cell],
                                     quad_f,
                                     timeStepper,
                                     par,
                                     problem_f[i_cell],
                                     solution,
                                     i_cell,
                                     mesh_collection.mesh.n_cell)

        # Compute time derivative of basis functions via finite differencing
        FiniteDiff.central!(solution.phi_1_t[i_cell], solution.phi_1[i_cell], par.dt)
        FiniteDiff.central!(solution.phi_2_t[i_cell], solution.phi_2[i_cell], par.dt)
        FiniteDiff.central!(solution.phi_3_t[i_cell], solution.phi_3[i_cell], par.dt)
    end
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    

    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    # Set up time integrator
    timeStepper = TimeIntegrator.ImplEuler(dof_collection.dof,
                                            mesh_collection.mesh,
                                            problem)
 
    solve_problem_MsFEM!(mesh_collection,
                           ref_el,
                           dof_collection.dof,
                           quad_f,
                           timeStepper,
                           par,
                           problem,
                           solution)
    # ---------------------------------------------------------------------
    # ---------------------------------------------------------------------
    
    return solution, mesh_collection
end


# -------------------------------------------------------------------------------------------


function solve_problem_MsFEM_basis!(mesh :: Mesh.TriangleMesh.TriMesh,
                                    ref_el :: FEM.AbstractRefEl,
                                    dof :: FEM.AbstractDof,
                                    quad :: Quad.AbstractQuad,
                                    timeStepper :: TimeIntegrator.AbstractTimeIntegrator,
                                    par :: Parameter.Parameter_MsFEM,
                                    problem :: Problem.AbstractBasisProblem,
                                    solution :: FEM.AbstractSolution,
                                    ind_cell :: Int64, n_cell :: Int64)

    println("\n\n--------------------------------------------------------------")
    println("Computing local multiscale FEM solution:   $(ind_cell) / $(n_cell)\n")
    println("\t Mesh type:   $(mesh.mesh_info)")
    if dof.is_periodic
        println("\t Problem type:   periodic --- $(problem.type_info)")
    else
        println("\t Problem type:   non-periodic --- $(problem.type_info)")
    end
    println("\t Time intergrator:   $(typeof(timeStepper))")
    println("\t number of time steps:   $(par.n_steps)")
    println("\t Number of active dofs:   $(dof.n_true_dof - dof.n_node_dirichlet)")
    println("\t Number of constraint dofs:   $(dof.n_node_dirichlet)")
    
    
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    
    # Setup initial data
    u_tmp = zeros(mesh.n_point, 3, par.n_steps+1)
    u_tmp[:,:,1] = Problem.u_init(problem, mesh.point)
    solution.phi_1[ind_cell][:,1] = u_tmp[:,1,1]
    solution.phi_2[ind_cell][:,1] = u_tmp[:,2,1]
    solution.phi_3[ind_cell][:,1] = u_tmp[:,3,1]

    # Make step from k_time to k_time+1
    M, A, D = [], [], []
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
          solution.mass[ind_cell] = M
        end
        
        if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
          A = FEM.assemble_advection(mesh,
                                       dof,
                                       ref_el,
                                       quad,
                                       problem,
                                       Time+par.dt)
          solution.advection[ind_cell,k_time+1] = A
        end

        if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
          D = FEM.assemble_diffusion(mesh,
                                       dof,
                                       ref_el,
                                       quad,
                                       problem,
                                       Time+par.dt)
          solution.diffusion[ind_cell,k_time+1] = D
        end

        # Zero forcing
        f = zeros(mesh.n_point,3)

        # Set the system matrices 
        TimeIntegrator.updateSystem!(timeStepper.systemData, 
                                      M, D-A, f, 
                                      dof, 
                                      u_tmp[:,:,k_time], 
                                      par.dt)

        # Make a single time step        
        TimeIntegrator.makeStep!(timeStepper, dof, view(u_tmp, :, :, k_time+1), u_tmp[:,:,k_time])
        solution.phi_1[ind_cell][:,k_time+1] = u_tmp[:,1,k_time+1]
        solution.phi_2[ind_cell][:,k_time+1] = u_tmp[:,2,k_time+1]
        solution.phi_3[ind_cell][:,k_time+1] = u_tmp[:,3,k_time+1]

        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")

    return nothing
end


# -------------------------------------------------------------------------------------------


function solve_problem_MsFEM!(mesh_collection :: Mesh.TriMesh_collection,
                              ref_el :: FEM.AbstractRefEl,
                              dof :: FEM.AbstractDof,
                              quad :: Quad.AbstractQuad,
                              timeStepper :: TimeIntegrator.AbstractTimeIntegrator,
                              par :: Parameter.Parameter_MsFEM,
                              problem :: Problem.AbstractProblem,
                              solution :: FEM.Solution_MsFEM)

    println("\n\n--------------------------------------------------------------")

    println("Computing multiscale FEM solution:\n")
    
    println("\t Mesh type:   $(mesh_collection.mesh.mesh_info)")
    if dof.is_periodic
        println("\t Problem type:   periodic --- $(problem.type_info)")
    else
        println("\t Problem type:   non-periodic --- $(problem.type_info)")
    end
    println("\t Time intergrator:   $(typeof(timeStepper))")
    println("\t number of coarse elements:   $(mesh_collection.mesh.n_cell)")
    println("\t average number of fine of elements:   $(mean(mesh_collection.n_elem_f))")
    println("\t number of all (fine) elements:   $(sum(mesh_collection.n_elem_f))")
    println("\t number of time steps:   $(par.n_steps)")
    println("\t Number of active dofs:   $(dof.n_true_dof - dof.n_node_dirichlet)")
    println("\t Number of constraint dofs:   $(dof.n_node_dirichlet)")
    
    # Setup initial data
    solution.u[:,1] = Problem.u_init(problem, mesh_collection.mesh.point)

    # Make step from k_time to k_time+1
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    for k_time in 1:par.n_steps
        # Current time index
        Time = k_time*par.dt

        if problem.is_transient_velocity
            k_adv = k_time + 1
        else
            k_adv = 2
        end

        if problem.is_transient_diffusion
            k_diff = k_time + 1
        else
            k_diff = 2
        end

        M = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_mass!(M,
                           solution.mass,
                           solution,
                           mesh_collection.mesh,
                           dof,
                           ref_el,
                           problem,
                           k_time+1)

        A = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_advection!(A,
                                solution.advection[:,k_adv],
                                solution,
                                mesh_collection.mesh,
                                dof,
                                ref_el,
                                problem,
                                k_time+1)

        # This line simply adds the derivative of the mass matrix term to
        # the advection matrix. This is done in place.
        FEM.assemble_mass_t!(A,
                             solution.mass,
                             solution,
                             mesh_collection.mesh,
                             dof,
                             ref_el,
                             problem,
                             k_time+1)

        D = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_diffusion!(D,
                                solution.diffusion[:,k_diff],
                                solution,
                                mesh_collection.mesh,
                                dof,
                                ref_el,
                                problem,
                                k_time+1)
        # Zero forcing
        f = zeros(mesh_collection.mesh.n_point)

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


