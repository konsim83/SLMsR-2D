function updateSystem!(systemData :: ImplEulerData,
                        M :: SparseMatrixCSC{Float64,Int64},
                        A :: SparseMatrixCSC{Float64,Int64},
                        f :: Array{Float64},
                        uOld :: Array{Float64},
                        dt :: Float64)
    
    # Just an abbreviation
    innd = systemData.ind_node_non_dirichlet
    ind = systemData.ind_node_dirichlet

    uOldDof = FEM.map_vec_mesh2dof(dof, uOld[:,:])
    fDof = FEM.map_vec_mesh2dof(dof, f[:,:])

    systemData.system_matrix[:,:] = (M - dt*A)[innd,innd]
    
    warn("Here is a little inconsistency.")
    
    systemData.system_rhs[:,:] = M[innd,innd] * uOldDof[innd,:] + fDof[innd]
    
    if !isempty(ind)
        systemData.system_rhs[:,:] +=  M[innd,ind] * uOldDof[ind,:] - (M - dt*A)[innd,ind] * uOldDof[ind,:]
    end
    

end






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

    Update the system at time index (k_time+1)*par.dt because we need
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
    system_data.mass[:,:] = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    FEM.assemble_mass!(system_data.mass,
                       solution.mass,
                       solution,
                       mesh,
                       dof,
                       ref_el,
                       par,
                       problem,
                       (k_time+1)*par.dt)

    system_data.advection[:,:] = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    FEM.assemble_advection!(system_data.advection,
                            solution.advection[:,k_adv],
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            par,
                            problem,
                            (k_time+1)*par.dt)

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
                         (k_time+1)*par.dt)

    system_data.diffusion[:,:] = sparse(dof.ind_test, dof.ind, zeros(Float64, length(dof.ind)), dof.n_true_dof, dof.n_true_dof)
    FEM.assemble_diffusion!(system_data.diffusion,
                            solution.diffusion[:,k_diff],
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            par,
                            problem,
                            (k_time+1)*par.dt)
    # ----   This is the fast version   ---
    
    if dof.n_node_neumann > 0
        error("Neumann boundary integral not implemented yet.")
    end

    system_data.system_matrix[:,:] = ( (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[system_data.ind_node_non_dirichlet,system_data.ind_node_non_dirichlet] )
    
    if dof.is_periodic
        system_data.system_rhs[:] = system_data.mass[system_data.ind_node_non_dirichlet,system_data.ind_node_non_dirichlet] * (FEM.map_vec_mesh2dof(dof, solution.u[:,k_time])[system_data.ind_node_non_dirichlet])
    else
        
        system_data.system_rhs[:] = (   (system_data.mass[system_data.ind_node_non_dirichlet,:]
                                      * solution.u[:,k_time] -

                                      (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[system_data.ind_node_non_dirichlet,system_data.ind_node_dirichlet]
                                      * solution.u[system_data.ind_node_dirichlet,k_time]) )
    end
    
end
# ----------------------------------
