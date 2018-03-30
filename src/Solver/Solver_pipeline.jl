function solve_problem!(mesh :: Mesh.TriangleMesh.TriMesh,
                        ref_el :: FEM.AbstractRefEl,
                        dof :: FEM.AbstractDof,
                        quad :: Quad.AbstractQuad,
                        timeStepper :: TimeIntegrator.AbstractTimeIntegrator,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractProblem,
                        solution :: FEM.Solution_FEM)

    println("\n\n--------------------------------------------------------------")
    println("Computing standard FEM solution:\n")
    
    println("\t Mesh type:   $(mesh.mesh_info)")
    if dof.is_periodic
        println("\t Problem type:   periodic --- $(problem.type_info)")
    else
        println("\t Problem type:   non-periodic --- $(problem.type_info)")
    end
    println("\t Time intergrator:   $(typeof(timeStepper))")
    println("\t number of elements:   $(mesh.n_cell)")
    println("\t number of time steps:   $(par.n_steps)")
    println("\t Number of active dofs:   $(dof.n_true_dof - dof.n_node_dirichlet)")
    println("\t Number of constraint dofs:   $(dof.n_node_dirichlet)")
    
    # Setup initial data
    solution.u[:,1] = Problem.u_init(problem, mesh.point)

    # Make step from k_time to k_time+1
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    for k_time in 1:par.n_steps
        # Current time index
        Time = k_time*par.dt


        # Assemble the system matrices only when necessary
        if k_time==1
          # Zero forcing
          f = zeros(size(D,1))

          M = FEM.assemble_mass(mesh,
                                 dof,
                                 ref_el,
                                 quad,
                                 problem)
        end
        
        if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
          A = FEM.assemble_advection(mesh,
                                       dof,
                                       ref_el,
                                       quad,
                                       par,
                                       problem,
                                       Time+par.dt)
        end

        if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
          D = FEM.assemble_diffusion(mesh,
                                       dof,
                                       ref_el,
                                       quad,
                                       par,
                                       problem,
                                       Time+par.dt)
        end

        

        # Set the system matrices 
        TimeIntegrator.updateSystem(timeStepper.systemData, M, D-A, f, par.dt)

        # Make a single time step
        TimeIntegrator.makeStep!(timeStepper, view(solution.u, :, k_time+1), solution.u[:,k_time], par.dt)

        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")

    return nothing
end


# -------------------------------------------------------------------------------------------


function solve_problem_local!(mesh :: Mesh.TriangleMesh.TriMesh,
                              ref_el :: FEM.AbstractRefEl,
                              dof :: FEM.AbstractDof,
                              quad :: Quad.AbstractQuad,
                              time_stepper :: TimeIntegrator.AbstractTimeIntegrator,
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
    N = par.n_steps
    p = Progress(N, 0.01, "Progress of time stepping...", 10)
    for k_time in 1:par.n_steps
        # Current time index
        Time = k_time*par.dt


        # Assemble the system matrices only when necessary
        if k_time==1
          # Zero forcing
          f = zeros(size(D,1),3)

          solution.mass[ind_cell] = FEM.assemble_mass(mesh,
                                             dof,
                                             ref_el,
                                             quad,
                                             problem)
        end
        
        if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
          solution.advection[ind_cell,k_time+1] = FEM.assemble_advection(mesh,
                                                       dof,
                                                       ref_el,
                                                       quad,
                                                       par,
                                                       problem,
                                                       Time+par.dt)
        end

        if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
          solution.diffusion[ind_cell,k_time+1] = FEM.assemble_diffusion(mesh,
                                                       dof,
                                                       ref_el,
                                                       quad,
                                                       par,
                                                       problem,
                                                       Time+par.dt)
        end

        

        # Set the system matrices 
        TimeIntegrator.updateSystem(timeStepper.systemData, M, D-A, f, par.dt)

        # Make a single time step
        error("How to pass the basis?")
        TimeIntegrator.makeStep!(timeStepper, view(solution.u, :, k_time+1), solution.u[:,k_time], par.dt)

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
                        time_stepper :: TimeIntegrator.AbstractTimeIntegrator,
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
    println("\t Time intergrator:   $(typeof(time_stepper))")
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
    
        error("Change to explicit non-inplace assembly.")

        M = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_mass!(M,
                           solution.mass,
                           solution,
                           mesh,
                           dof,
                           ref_el,
                           par,
                           problem,
                           k_time+1)

        A = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_advection!(A,
                                solution.advection[:,k_adv],
                                solution,
                                mesh,
                                dof,
                                ref_el,
                                par,
                                problem,
                                k_time+1)

        # This line simply adds the derivative of the mass matrix term to
        # the advection matrix. This is done in place.
        FEM.assemble_mass_t!(A,
                             solution.mass,
                             solution,
                             mesh,
                             dof,
                             ref_el,
                             par,
                             problem,
                             k_time+1)

        D = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_diffusion!(D,
                                solution.diffusion[:,k_diff],
                                solution,
                                mesh,
                                dof,
                                ref_el,
                                par,
                                problem,
                                k_time+1)
        # Zero forcing
        f = zeros(size(D,1))

        # Set the system matrices 
        TimeIntegrator.updateSystem(timeStepper.systemData, M, D-A, f, par.dt)

        # Make a single time step
        TimeIntegrator.makeStep!(timeStepper, view(solution.u, :, k_time+1), solution.u[:,k_time], par.dt)

        next!(p)
    end # end for
    println("..... done.")
    println("\n--------------------------------------------------------------\n\n")

    return nothing
end


# -------------------------------------------------------------------------------------------