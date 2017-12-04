function solve_problem!(mesh :: Mesh.TriMesh,
                        ref_el :: FEM.AbstractRefEl,
                        dof :: FEM.AbstractDof,
                        quad :: Quad.AbstractQuad,
                        time_stepper :: Time_integrator.AbstractTime_integrator,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractProblem,
                        solution :: FEM.AbstractSolution)

    """

    This is the actual solution routine. The pipeline is

    1. Set up initial data

    2. Set up system that can make a step from k to k+1

    3. Call the actual time stepper. Note that the type system_data
    must be created accoring to the type of time_stepper (it is
    contained in it as a subtype).

    """

    # Setup initial data
    solution.u[:,1] = Problem.u_init(problem, mesh.point)

    # Make step from k_time to k_time+1
    for k_time=1:par.n_steps        
        Time_integrator.make_step!(solution,
                                   time_stepper,
                                   mesh,
                                   dof,
                                   ref_el,
                                   quad,
                                   par,
                                   problem,
                                   k_time)
        display("--------------------------------------")
        display(solution.u[:,k_time+1])
    end # end for    
end
