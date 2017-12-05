function solve_FEM(par :: Parameter.Parameter_FEM, problem :: T) where {T<:Problem.AbstractProblem}

    # Build mesh of unit square (0,1)x(0,1)
    mesh = Mesh.mesh_unit_square(par.n_edge_per_seg)

    # Set up reference element
    ref_el = FEM.RefEl_Pk{par.n_order_FEM}()

    # Set up quadrature rule
    quad = Quad.Quad_simplex(3)

    # Set up degrees of freedom handler
    dof = FEM.Dof_Pk_periodic_square{ref_el.n_order}(mesh);

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
