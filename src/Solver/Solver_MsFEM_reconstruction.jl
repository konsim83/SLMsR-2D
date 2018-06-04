function solve_MsFEM_periodic_square_reconstruction(par :: Parameter.Parameter_MsFEM,
                                                    problem :: T) where {T<:Problem.AbstractProblem}
    
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
    problem_f = Array{Problem.AbstractBasisProblem,1}(mesh_collection.mesh.n_cell)
    for i_cell in 1:mesh_collection.mesh.n_cell
        # Set up local problem by geometry
        tri = Geometry.Triangle(FEM.map_ref_point(dof, ref_el.node, i_cell))
        problem_f[i_cell] = Problem.BasisFun(problem, tri)
    end
    dof_collection = FEM.Dof_collection(mesh_collection, dof, problem_f, ref_el_f.n_order)

    # Set up solution structure
    solution = FEM.Solution_MsFEM(dof_collection, par)


    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    # Setup initial data for basis functions
    timeStepper_local = Array{TimeIntegrator.ImplEuler,1}(mesh_collection.mesh.n_cell)
    for ind_cell in 1:mesh_collection.mesh.n_cell
        u_init_tmp = Problem.u_init(problem_f[ind_cell], mesh_collection.mesh_f[ind_cell].point)
        solution.phi_1[ind_cell][:,1] = u_init_tmp[:,1]
        solution.phi_2[ind_cell][:,1] = u_init_tmp[:,2]
        solution.phi_3[ind_cell][:,1] = u_init_tmp[:,3]

        timeStepper_local[ind_cell] = TimeIntegrator.ImplEuler(dof_collection.dof_f[ind_cell],
                                                                    mesh_collection.mesh_f[ind_cell],
                                                                    problem_f[ind_cell])
    end

    # Setup initial data
    solution.u[1:mesh_collection.mesh.n_point,1] = Problem.u_init(problem, mesh_collection.mesh.point)
    
    timeStepper = TimeIntegrator.ImplEuler(dof_collection.dof,
                                            mesh_collection.mesh,
                                            problem)


    println("\n\n--------------------------------------------------------------")

    println("Computing reconstruction multiscale FEM solution:\n")
    
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


    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    for k_time=1:par.n_steps
        # Time at index k
        Time = (k_time-1)*par.dt

        # Reconstruct the next basis function (due to implicit scheme)
        # Reconstruction.SemiLagrange_L2_opt!(solution, 
        #                                         timeStepper_local,
        #                                         mesh_collection,
        #                                         dof_collection,
        #                                         ref_el_f,
        #                                         quad_f,
        #                                         par,
        #                                         problem,
        #                                         problem_f,
        #                                         Time + par.dt, k_time + 1)

        Reconstruction.SemiLagrange_L2_opt_new!(solution, 
                                                timeStepper_local,
                                                mesh_collection,
                                                dof_collection,
                                                ref_el_f,
                                                quad_f,
                                                par,
                                                problem,
                                                problem_f,
                                                Time + par.dt, k_time + 1)

        M, Mt = FEM.assemble_mass(solution,
                                   mesh_collection,
                                   dof_collection,
                                   ref_el, ref_el_f,
                                   quad_f,
                                   problem,
                                   k_time+1)

        A = FEM.assemble_advection(solution,
                                   mesh_collection,
                                   dof_collection,
                                   ref_el, ref_el_f,
                                   quad_f,
                                   problem,
                                   k_time+1, Time + par.dt)

        D = FEM.assemble_diffusion(solution,
                                   mesh_collection,
                                   dof_collection,
                                   ref_el, ref_el_f,
                                   quad_f,
                                   problem,
                                   k_time+1, Time + par.dt)

        # Zero forcing
        f = zeros(mesh_collection.mesh.n_point)

        # Set the system matrices 
        TimeIntegrator.updateSystem!(timeStepper.systemData, M, D-A-Mt, f, dof, solution.u[:,k_time], par.dt)

        # Make a single time step
        TimeIntegrator.makeStep!(timeStepper, dof, view(solution.u, :, k_time+1), solution.u[:,k_time])
        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")
    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    
    return solution, mesh_collection
end