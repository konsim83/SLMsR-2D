function solve_FEM_periodic_square(par :: Parameter.Parameter_FEM, problem :: T) where {T<:Problem.AbstractProblem}

    # Build mesh of unit square (0,1)x(0,1)
    mesh = Mesh.mesh_unit_square(par.n_edge_per_seg)

    # Set up reference element
    ref_el = FEM.RefEl_Pk{par.n_order_FEM}()

    # Set up quadrature rule
    quad = Quad.Quad_simplex(par.n_order_quad)

    # Set up degrees of freedom handler
    dof = FEM.Dof_Pk_periodic_square{ref_el.n_order}(mesh)

    # Set up solution structure
    solution = FEM.Solution_FEM(dof, par)

    # Set up time integrator
    time_stepper = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(dof, mesh, problem)
    # +++++++++++++++++++

    # Call actual solver. Pass solution data by reference.       
    solve_problem!(mesh,
                   ref_el,
                   dof,
                   quad,
                   time_stepper,
                   par,
                   problem,
                   solution)
    
    return solution, mesh
end


# ------------------------------------------------------------------------------------------------------------------------------



function solve_FEM_simplex(par :: Parameter.Parameter_FEM, problem :: T) where {T<:Problem.AbstractProblem}

    # Build mesh of unit square (0,1)x(0,1)
    mesh = Mesh.mesh_unit_simplex_uniform_edges(par.n_edge_per_seg)

    # Set up reference element
    ref_el = FEM.RefEl_Pk{par.n_order_FEM}()

    # Set up quadrature rule
    quad = Quad.Quad_simplex(par.n_order_quad)

    # Set up degrees of freedom handler
    dof = FEM.Dof_Pk{ref_el.n_order}(mesh);

    # Set up solution structure
    solution = FEM.Solution_FEM(dof, par)

    # Set up time integrator
    time_stepper = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(dof, mesh, problem)
    # +++++++++++++++++++

    # Call actual solver. Pass solution data by reference.       
    solve_problem!(mesh,
                   ref_el,
                   dof,
                   quad,
                   time_stepper,
                   par,
                   problem,
                   solution)
    
    return solution, mesh
end



function solve_MsFEM_periodic_square(par :: Parameter.Parameter_MsFEM, problem :: T) where {T<:Problem.AbstractProblem}

    # Build mesh of unit square (0,1)x(0,1)
    m_coarse = Mesh.mesh_unit_square(par.n_edge_per_seg)
    mesh_collection = Mesh.TriMesh_collection(m_coarse, par.n_edge_per_seg_f)
    #=
    mesh_simplex = Mesh.mesh_unit_simplex_uniform_edges(par.n_edge_per_seg)
    mesh_collection = Mesh.Trimesh_collection(m_coarse, mesh_simplex)
    =#
    

    # Set up reference element
    ref_el = FEM.RefEl_Pk{1}()
    ref_el_f = FEM.RefEl_Pk{par.n_order_FEM_f}()

    # Set up quadrature rule
    quad_f = Quad.Quad_simplex(par.n_order_quad_f)

    # Set up degrees of freedom handler
    dof_collection = FEM.Dof_collection{ref_el_f.n_order}(mesh_collection)

    # Set up solution structure
    solution = FEM.Solution_MsFEM(dof_collection, par)


    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # Solve the local problems
    for i_cell in 1:mesh_collection.mesh.n_cell
        # Set up local problem by geometry
        tri = Geometry.Triangle(Mesh.map_ref_point(mesh_collection.mesh, ref_el.node, i_cell))
        problem_f = Problem.BasisFun(problem, tri)
        
        # Set up time integrator
        time_stepper = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(dof_collection.dof_f[i_cell],
                                                                                            mesh_collection.mesh_f[i_cell],
                                                                                            problem_f)
        # Call actual solver. Pass solution data by reference.  Solves
        # the i-th cell problem.
        solve_problem_local!(mesh_collection.mesh_f[i_cell],
                             ref_el_f,
                             dof_collection.dof_f[i_cell],
                             quad_f,
                             time_stepper,
                             par,
                             problem_f,
                             solution,
                             i_cell)
    end
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    

    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    #=
    # Set up time integrator
    time_stepper = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(dof_collection.dof,
                                                                                        mesh_collection.mesh,
                                                                                        problem)

    # Call actual solver. Pass solution data by reference.       
    solve_problem!(mesh_collection.mesh,
                   ref_el,
                   dof_collection.dof,
                   time_stepper,
                   par,
                   problem,
                   solution)
    =#
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    return solution, mesh_collection
end



