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
    solve_problem!(mesh,
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
    problem_f = Array{Problem.AbstractBasisProblem,1}(mesh_collection.mesh.n_cell)
    for i_cell in 1:mesh_collection.mesh.n_cell
        # Set up local problem by geometry
        tri = Geometry.Triangle(FEM.map_ref_point(dof, ref_el.node, i_cell))
        problem_f[i_cell] = Problem.BasisFun(problem, tri)
    end
    dof_collection = FEM.Dof_collection(mesh_collection, dof, problem_f, ref_el_f.n_order)

    # Set up solution structure
    solution = FEM.Solution_MsFEM(dof_collection, par)


    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # Solve the local problems
    for i_cell in 1:mesh_collection.mesh.n_cell        
        # Set up time integrator
        timeStepper = TimeIntegrator.ImplEuler(dof_collection.dof_f[i_cell],
                                                mesh_collection.mesh_f[i_cell],
                                                problem_f[i_cell])
        # Call actual solver. Pass solution data by reference.  Solves
        # the i-th cell problem.
        solve_problem_local!(mesh_collection.mesh_f[i_cell],
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
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # Set up time integrator
    timeStepper = TimeIntegrator.ImplEuler{TimeIntegrator.System_data_implEuler_ADE}(dof_collection.dof,
                                                                                        mesh_collection.mesh,
                                                                                        problem)
 
    solve_problem!(mesh_collection,
                   ref_el,
                   dof_collection.dof,
                   quad_f,
                   timeStepper,
                   par,
                   problem,
                   solution)
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    return solution, mesh_collection
end


# ------------------------------------------------------------------------------------------------------------------------------


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
    for ind_cell in 1:mesh_collection.mesh.n_cell
        u_init_tmp = Problem.u_init(problem_f[ind_cell], mesh_collection.mesh_f[ind_cell].point)
        solution.phi_1[ind_cell][:,1] = u_init_tmp[:,1]
        solution.phi_2[ind_cell][:,1] = u_init_tmp[:,2]
        solution.phi_3[ind_cell][:,1] = u_init_tmp[:,3]
    end

    # Setup initial data
    solution.u[1:mesh_collection.mesh.n_point,1] = Problem.u_init(problem, mesh_collection.mesh.point)

    # Reconstruct the basis at time step k_time+1
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    for k_time=1:par.n_steps
        # Time at index k
        Time = (k_time-1)*par.dt

        Reconstruction.SemiLagrange_L2_opt(solution, 
                                            mesh_collection,
                                            dof_collection,
                                            par,
                                            problem,
                                            problem_f,
                                            Time + par.dt, k_time + 1)
        error("Now integrate in time...")
        TimeIntegrator.make_step!(solution,
                                   timeStepper,
                                   mesh,
                                   dof,
                                   ref_el,
                                   quad,
                                   par,
                                   problem,
                                   k_time,
                                   ind_cell)
        next!(p)
    end # end for
    # ----------------------------------------------------------------
    # ----------------------------------------------------------------
    
    return solution, mesh_collection
end