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

    println("\n\n--------------------------------------------------------------")
    println("Computing standard FEM solution:\n")
    println("\t Mesh type:   $(mesh.mesh_info)")
    if dof.is_periodic
        println("\t Problem type:   periodic --- $(problem.type_info)")
    else
        println("\t Problem type:   non-periodic --- $(problem.type_info)")
    end
    println("\t Time intergrator:   $(typeof(time_stepper))")
    println("\t number of time steps:   $(par.n_steps)")
    println("\t Number of active dofs:   $(dof.n_true_dof - dof.n_node_dirichlet)")
    println("\t Number of constraint dofs:   $(dof.n_node_dirichlet)")
    
    
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    
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
        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")

    return nothing
end
