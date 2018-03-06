# ----------------------------------
function update_system!(solution :: FEM.AbstractSolution,
                        system_data :: System_data_implEuler_ADE,
                        mesh :: Mesh.TriangleMesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.AbstractRefEl,
                        quad :: Quad.AbstractQuad,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractPhysicalProblem,
                        k_time :: Int64)

    """

    Update the system at time index k_time+1 because we need
    information at the next time step.

    """    
    
    # ----   This is the fast version   ---
    # if k_time==1
    #     system_data.mass = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    #     FEM.assemble_mass!(system_data.mass,
    #                        mesh,
    #                        dof,
    #                        ref_el,
    #                        quad,
    #                        par,
    #                        problem)
    # end

    # if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
    #     system_data.advection = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    #     FEM.assemble_advection!(system_data.advection,
    #                             mesh,
    #                             dof,
    #                             ref_el,
    #                             quad,
    #                             par,
    #                             problem,
    #                             k_time+1)
    # end

    # if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
    #     system_data.diffusion = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    #     FEM.assemble_diffusion!(system_data.diffusion,
    #                             mesh,
    #                             dof,
    #                             ref_el,
    #                             quad,
    #                             par,
    #                             problem,
    #                             k_time+1)
    # end
    # ----   This is the fast version   ---
    


    
    # ----   This is the slow version   ---
    if k_time==1
    system_data.mass = FEM.assemble_mass(mesh,
                                         dof,
                                         ref_el,
                                         quad,
                                         par,
                                         problem)
    end
    
    if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
    system_data.advection = FEM.assemble_advection(mesh,
                                                   dof,
                                                   ref_el,
                                                   quad,
                                                   par,
                                                   problem,
                                                   k_time+1)
    end

    if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
    system_data.diffusion = FEM.assemble_diffusion(mesh,
                                                   dof,
                                                   ref_el,
                                                   quad,
                                                   par,
                                                   problem,
                                                   k_time+1)
    end
    # ----   This is the slow version   ---
    

    
    if dof.n_node_neumann > 0
        error("Neumann boundary integral not implemented yet.")
    end

    system_data.system_matrix = ( (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] )
    
    if dof.is_periodic
        system_data.system_rhs = system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] * (FEM.map_vec_mesh2dof(dof, solution.u[:,k_time])[dof.ind_node_non_dirichlet])
    else
        
        system_data.system_rhs = (   (system_data.mass[dof.ind_node_non_dirichlet,:]
                                      * solution.u[:,k_time] -

                                      (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_dirichlet]
                                      * solution.u[dof.ind_node_dirichlet,k_time]) )
    end
    
end
# ----------------------------------


# ----------------------------------
function update_system!(solution :: FEM.AbstractSolution,
                        system_data :: System_data_implEuler_ADE,
                        mesh :: Mesh.TriangleMesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.AbstractRefEl,
                        quad :: Quad.AbstractQuad,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractBasisProblem,
                        k_time :: Int64,
                        ind_cell :: Int64)

    """

    Update the system at time index k_time+1 because we need
    information at the next time step.

    """
    
    # ----   This is the fast version   ---
    if k_time==1
        system_data.mass = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_mass!(system_data.mass,
                           mesh,
                           dof,
                           ref_el,
                           quad,
                           par,
                           problem)
        solution.mass[ind_cell][:,:] = copy(system_data.mass)
    end

    if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
        system_data.advection = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_advection!(system_data.advection,
                                mesh,
                                dof,
                                ref_el,
                                quad,
                                par,
                                problem,
                                k_time+1)
        solution.advection[ind_cell, k_time+1][:,:] = copy(system_data.advection)
    end

    if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
        system_data.diffusion = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
        FEM.assemble_diffusion!(system_data.diffusion,
                                mesh,
                                dof,
                                ref_el,
                                quad,
                                par,
                                problem,
                                k_time+1)
        solution.diffusion[ind_cell, k_time+1][:,:] = copy(system_data.diffusion)
    end
    # ----   This is the fast version   ---
    


    #=
    # ----   This is the slow version   ---
    if k_time==1
    system_data.mass = FEM.assemble_mass(mesh,
                                         dof,
                                         ref_el,
                                         quad,
                                         par,
                                         problem)
    end
    
    if (k_time==1) || (k_time>=1 && problem.is_transient_velocity)
    system_data.advection = FEM.assemble_advection(mesh,
                                                   dof,
                                                   ref_el,
                                                   quad,
                                                   par,
                                                   problem,
                                                   k_time+1)
    end

    if (k_time==1) || (k_time>=1 && problem.is_transient_diffusion)
    system_data.diffusion = FEM.assemble_diffusion(mesh,
                                                   dof,
                                                   ref_el,
                                                   quad,
                                                   par,
                                                   problem,
                                                   k_time+1)
    end
    # ----   This is the slow version   ---
    =#

    
    if dof.n_node_neumann > 0
        error("Neumann boundary integral not implemented yet.")
    end

    system_data.system_matrix = ( (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] )
    
    if dof.is_periodic
        system_data.system_rhs = system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] *
            (FEM.map_ind_mesh2dof(dof, solution.phi[ind_cell][:,k_time])[dof.ind_node_non_dirichlet])
    else

        phi_tmp =  [solution.phi_1[ind_cell][:,k_time] solution.phi_2[ind_cell][:,k_time] solution.phi_3[ind_cell][:,k_time]]
        system_data.system_rhs = (   (system_data.mass[dof.ind_node_non_dirichlet,:]
                                      * phi_tmp -

                                      (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_dirichlet]
                                      * phi_tmp[dof.ind_node_dirichlet,:]) )
    end
    
end
# ----------------------------------


# ----------------------------------
function update_system!(solution :: FEM.Solution_MsFEM,
                        system_data :: System_data_implEuler_ADE,
                        mesh :: Mesh.TriangleMesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.AbstractRefEl,
                        quad :: Quad.AbstractQuad,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractPhysicalProblem,
                        k_time :: Int64)

    """

    Update the system at time index k_time+1 because we need
    information at the next time step.

    """

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
    
    # ----   This is the fast version   ---
    system_data.mass = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    FEM.assemble_mass!(system_data.mass,
                       solution.mass,
                       solution,
                       mesh,
                       dof,
                       ref_el,
                       par,
                       problem,
                       k_time+1)

    system_data.advection = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    FEM.assemble_advection!(system_data.advection,
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
    FEM.assemble_mass_t!(system_data.advection,
                         solution.mass,
                         solution,
                         mesh,
                         dof,
                         ref_el,
                         par,
                         problem,
                         k_time+1)

    system_data.diffusion = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    FEM.assemble_diffusion!(system_data.diffusion,
                            solution.diffusion[:,k_diff],
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            par,
                            problem,
                            k_time+1)
    # ----   This is the fast version   ---
    
    if dof.n_node_neumann > 0
        error("Neumann boundary integral not implemented yet.")
    end

    system_data.system_matrix = ( (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] )
    
    if dof.is_periodic
        system_data.system_rhs = system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] * (FEM.map_ind_mesh2dof(dof, solution.u[:,k_time])[dof.ind_node_non_dirichlet])
    else
        
        system_data.system_rhs = (   (system_data.mass[dof.ind_node_non_dirichlet,:]
                                      * solution.u[:,k_time] -

                                      (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_dirichlet]
                                      * solution.u[dof.ind_node_dirichlet,k_time]) )
    end
    
end
# ----------------------------------
