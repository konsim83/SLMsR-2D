function solve_problem!(mesh :: Mesh.TriMesh,
                        ref_el :: FEM.RefEl,
                        dof :: FEM.Dof,
                        quad :: Quad.Quad_top,
                        par :: Parameter.Parameter_top,
                        problem :: Problem.Problem_top,
                        solution :: FEM.Solution)

    """

    This is the actual solution routine. The pipeline is

    1. Set up initial data

    2. Set up system that can make a step from k to k+1

    3. Call the actual time stepper. Note that the type setup_system
    must be created accoring to the type of time_stepper.

    """

    # Setup initial data
    solution.u[:,1] = Problem.u_init(problem, mesh.point)

    # Make step from k_time to k_time+1
    for k_time=1:par.n_steps
        system_data = Setup_system_ADE_implEuler(solution,
                                                 mesh,
                                                 dof,
                                                 ref_el,
                                                 quad,
                                                 par,
                                                 problem,
                                                 k_time)
        
        Time_integrator.make_step!(dof, system_data, solution)
    end # end for    
end