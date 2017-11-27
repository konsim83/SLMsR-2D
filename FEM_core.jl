function solve_FEM(par :: Parameter.Parameter_FEM, problem :: Problem.Problem_top)

    # Build mesh of unit square (0,1)x(0,1)
    mesh = Mesh.mesh_unit_square(par.n_edge_per_seg)

    # Set up reference element
    ref_el = RefEl_Pk{par.n_order_FEM}()

    # Set up quadrature rule
    quad = Quad.Quad_simplex(3)

    # Set up degrees of freedom handler
    dof = Dof_Pk_periodic_square{ref_el.n_order}(mesh);

    # Set up solution structure
    solution = Solution_FEM(dof, par)

    # Set up time integrator
    #time_stepper = Time_integrator.Implicit_Euler(par)
    # +++++++++++++++++++

    # Call actual solver. Pass solution data by reference.
    solve_problem!(mesh,
                   ref_el,
                   dof,
                   quad,
                   par,
                   problem,
                   solution)
    
    return solution, mesh
end
