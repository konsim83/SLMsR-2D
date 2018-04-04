# ---------------------------------------------------------------------------------------------
# Multiscale FEM assembling (direct)
# ---------------------------------------------------------------------------------------------

function assemble_mass(solution :: Solution_MsFEM,
                        mesh_collection :: Mesh.TriMesh_collection,
                        dof_collection :: FEM.Dof_collection,
                        ref_el :: FEM.RefEl_Pk{1},
                        ref_el_f :: FEM.RefEl_Pk{N},
                        quad_f :: Quad.AbstractQuad,
                        problem :: Problem.AbstractPhysicalProblem,
                        k_time :: Int64) where N
    

    m = Array{Float64,3}(ref_el.n_basis, ref_el.n_basis, mesh_collection.mesh.n_cell)
    mt = Array{Float64,3}(ref_el.n_basis, ref_el.n_basis, mesh_collection.mesh.n_cell)
    
    for k = 1:mesh_collection.mesh.n_cell
        Phi_t = [solution.phi_1_t[k][:,k_time] solution.phi_2_t[k][:,k_time] solution.phi_3_t[k][:,k_time]]
        Phi = [solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        mass = assemble_mass(mesh_collection.mesh_f[k],
                                        dof_collection.dof_f[k],
                                        ref_el_f,
                                        quad_f,
                                        problem)

        m[:,:,k] = Phi_test * mass * Phi
        mt[:,:,k] = Phi_test * mass * Phi_t
    end

    M = sparse(dof_collection.dof.ind_test, 
                        dof_collection.dof.ind, 
                        vec(m),
                        dof_collection.dof.n_true_dof,
                        dof_collection.dof.n_true_dof)

    Mt = sparse(dof_collection.dof.ind_test, 
                        dof_collection.dof.ind, 
                        vec(mt),
                        dof_collection.dof.n_true_dof,
                        dof_collection.dof.n_true_dof)
    
    return M, Mt
end



# **************************************************************************
# **************************************************************************
# **************************************************************************


function assemble_advection(solution :: Solution_MsFEM,
                                mesh_collection :: Mesh.TriMesh_collection,
                                dof_collection :: FEM.Dof_collection,
                                ref_el :: FEM.RefEl_Pk{1},
                                ref_el_f :: FEM.RefEl_Pk{N},
                                quad_f :: Quad.AbstractQuad,
                                problem :: Problem.AbstractPhysicalProblem,
                                k_time :: Int64, time_idx :: Float64) where N
    

    adv = Array{Float64,3}(ref_el.n_basis, ref_el.n_basis, mesh_collection.mesh.n_cell)
    
    for k = 1:mesh_collection.mesh.n_cell
        Phi = [solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        advection = assemble_advection(mesh_collection.mesh_f[k],
                                        dof_collection.dof_f[k],
                                        ref_el_f,
                                        quad_f,
                                        problem,
                                        time_idx)

        adv[:,:,k] = Phi_test * advection * Phi
    end

    Mat_global = sparse(dof_collection.dof.ind_test, 
                        dof_collection.dof.ind, 
                        vec(adv),
                        dof_collection.dof.n_true_dof,
                        dof_collection.dof.n_true_dof)
    
    return Mat_global
end



# **************************************************************************
# **************************************************************************
# **************************************************************************



function assemble_diffusion(solution :: Solution_MsFEM,
                                mesh_collection :: Mesh.TriMesh_collection,
                                dof_collection :: FEM.Dof_collection,
                                ref_el :: FEM.RefEl_Pk{1},
                                ref_el_f :: FEM.RefEl_Pk{N},
                                quad_f :: Quad.AbstractQuad,
                                problem :: Problem.AbstractPhysicalProblem,
                                k_time :: Int64, time_idx :: Float64) where N
    

    dif = Array{Float64,3}(ref_el.n_basis, ref_el.n_basis, mesh_collection.mesh.n_cell)
    
    for k = 1:mesh_collection.mesh.n_cell
        Phi = [solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        diffusion = assemble_diffusion(mesh_collection.mesh_f[k],
                                        dof_collection.dof_f[k],
                                        ref_el_f,
                                        quad_f,
                                        problem,
                                        time_idx)

        dif[:,:,k] = Phi_test * diffusion * Phi
    end

    Mat_global = sparse(dof_collection.dof.ind_test, 
                        dof_collection.dof.ind, 
                        vec(dif),
                        dof_collection.dof.n_true_dof,
                        dof_collection.dof.n_true_dof)
    
    return Mat_global
end