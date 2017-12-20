function solve_problem!(mesh :: Mesh.TriMesh,
                        ref_el :: FEM.AbstractRefEl,
                        dof :: FEM.AbstractDof,
                        quad :: Quad.AbstractQuad,
                        time_stepper :: Time_integrator.AbstractTime_integrator,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractProblem,
                        solution :: FEM.Solution_FEM)

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
    println("\t number of elements:   $(mesh.n_cell)")
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


# -------------------------------------------------------------------------------------------


function solve_problem_local!(mesh :: Mesh.TriMesh,
                              ref_el :: FEM.AbstractRefEl,
                              dof :: FEM.AbstractDof,
                              quad :: Quad.AbstractQuad,
                              time_stepper :: Time_integrator.AbstractTime_integrator,
                              par :: Parameter.Parameter_MsFEM,
                              problem :: Problem.AbstractBasisProblem,
                              solution :: FEM.AbstractSolution,
                              ind_cell :: Int64, n_cell :: Int64)

    """

    This is the actual solution routine. The pipeline is

    1. Set up initial data

    2. Set up system that can make a step from k to k+1

    3. Call the actual time stepper. Note that the type system_data
    must be created accoring to the type of time_stepper (it is
    contained in it as a subtype).

    """

    println("\n\n--------------------------------------------------------------")
    println("Computing local multiscale FEM solution:   $(ind_cell) / $(n_cell)\n")
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
    u_init_tmp = Problem.u_init(problem, mesh.point)
    solution.phi_1[ind_cell][:,1] = u_init_tmp[:,1]
    solution.phi_2[ind_cell][:,1] = u_init_tmp[:,2]
    solution.phi_3[ind_cell][:,1] = u_init_tmp[:,3]

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
                                   k_time,
                                   ind_cell)
        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")

    return nothing
end


# -------------------------------------------------------------------------------------------


function solve_problem!(mesh_collection :: Mesh.TriMesh_collection,
                        ref_el :: FEM.AbstractRefEl,
                        dof :: FEM.AbstractDof,
                        quad :: Quad.AbstractQuad,
                        time_stepper :: Time_integrator.AbstractTime_integrator,
                        par :: Parameter.Parameter_MsFEM,
                        problem :: Problem.AbstractProblem,
                        solution :: FEM.Solution_MsFEM)

    """

    This is the actual solution routine. The pipeline is

    1. Set up initial data

    2. Set up system that can make a step from k to k+1

    3. Call the actual time stepper. Note that the type system_data
    must be created accoring to the type of time_stepper (it is
    contained in it as a subtype).

    """

    println("\n\n--------------------------------------------------------------")

    println("Computing multiscale FEM solution:\n")
    
    println("\t Mesh type:   $(mesh_collection.mesh.mesh_info)")
    if dof.is_periodic
        println("\t Problem type:   periodic --- $(problem.type_info)")
    else
        println("\t Problem type:   non-periodic --- $(problem.type_info)")
    end
    println("\t Time intergrator:   $(typeof(time_stepper))")
    println("\t number of coarse elements:   $(mesh_collection.mesh.n_cell)")
    println("\t average number of fine of elements:   $(mean(mesh_collection.n_elem_f))")
    println("\t number of all (fine) elements:   $(sum(mesh_collection.n_elem_f))")
    println("\t number of time steps:   $(par.n_steps)")
    println("\t Number of active dofs:   $(dof.n_true_dof - dof.n_node_dirichlet)")
    println("\t Number of constraint dofs:   $(dof.n_node_dirichlet)")

    
    
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    
    # Setup initial data
    solution.u[:,1] = Problem.u_init(problem, mesh_collection.mesh.point)

    # Make step from k_time to k_time+1
    for k_time=1:par.n_steps
        Time_integrator.make_step!(solution,
                                   time_stepper,
                                   mesh_collection.mesh,
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


# -------------------------------------------------------------------------------------------
