# ---------------------------------------------------------------------------------------------
# Multiscale FEM assembling
# ---------------------------------------------------------------------------------------------


function assemble_mass!(M :: SparseMatrixCSC{Float64,Int64},
                        mass :: Array{SparseMatrixCSC{Float64,Int64},1},
                        solution :: Solution_MsFEM,
                        mesh :: Mesh.TriangleMesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.RefEl_Pk{1},
                        problem :: Problem.AbstractPhysicalProblem,
                        k_time :: Int64)

    m_loc = Array{Float64,2}(undef, ref_el.n_basis, ref_el.n_basis)
    
    for k = 1:mesh.n_cell
        Phi = [solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        m_loc[:,:] = Phi_test * mass[k] * Phi
        for l in 1:9
            M[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += m_loc[l]
        end
    end
    
    return nothing
end


# ---------------------------------------------------------------------------------------------


function assemble_mass_t!(Mt  :: SparseMatrixCSC{Float64,Int64},
                          mass :: Array{SparseMatrixCSC{Float64,Int64},1},
                          solution :: Solution_MsFEM,
                          mesh :: Mesh.TriangleMesh.TriMesh,
                          dof :: FEM.AbstractDof,
                          ref_el :: FEM.RefEl_Pk{1},
                          problem :: Problem.AbstractPhysicalProblem,
                          k_time :: Int64)

    mt_loc = Array{Float64,2}(undef, ref_el.n_basis, ref_el.n_basis)
    
    for k = 1:mesh.n_cell
        Phi_t = [solution.phi_1_t[k][:,k_time] solution.phi_2_t[k][:,k_time] solution.phi_3_t[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        mt_loc[:,:] = Phi_test * mass[k] * Phi_t
        for l in 1:9
            Mt[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += mt_loc[l]
        end
    end
    
    return nothing
end


# ---------------------------------------------------------------------------------------------


function assemble_advection!(A  :: SparseMatrixCSC{Float64,Int64},
                             advection :: Array{SparseMatrixCSC{Float64,Int64},1},
                             solution :: Solution_MsFEM,
                             mesh :: Mesh.TriangleMesh.TriMesh,
                             dof :: FEM.AbstractDof,
                             ref_el :: FEM.RefEl_Pk{1},
                             problem :: Problem.AbstractPhysicalProblem,
                             k_time :: Int64)

    a_loc = Array{Float64,2}(undef, ref_el.n_basis, ref_el.n_basis)
    
    for k = 1:mesh.n_cell
        Phi = [solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        a_loc[:,:] = Phi_test * advection[k] * Phi
        for l in 1:9
            A[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += a_loc[l]
        end
    end
    
    return nothing
end


# ---------------------------------------------------------------------------------------------


function assemble_diffusion!(D  :: SparseMatrixCSC{Float64,Int64},
                             diffusion :: Array{SparseMatrixCSC{Float64,Int64},1},
                             solution :: Solution_MsFEM,
                             mesh :: Mesh.TriangleMesh.TriMesh,
                             dof :: FEM.AbstractDof,
                             ref_el :: FEM.RefEl_Pk{1},
                             problem :: Problem.AbstractPhysicalProblem,
                             k_time :: Int64)
    
    d_loc = Array{Float64,2}(undef, ref_el.n_basis, ref_el.n_basis)
    
    for k = 1:mesh.n_cell
        Phi = [solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]]
        Phi_test = transpose([solution.phi_1[k][:,k_time] solution.phi_2[k][:,k_time] solution.phi_3[k][:,k_time]])
        
        d_loc[:,:] = Phi_test * diffusion[k] * Phi
        for l in 1:9
            D[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += d_loc[l]
        end
    end
    
    return nothing
end
