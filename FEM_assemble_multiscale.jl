# ---------------------------------------------------------------------------------------------
# Multiscale FEM assembling
# ---------------------------------------------------------------------------------------------


function assemble_mass!(M  :: SparseMatrixCSC{Float64,Int64},
                        solution :: Solution_MsFEM,
                        mesh :: Mesh.TriMesh,
                        dof :: FEM.Dof,
                        ref_el :: FEM.RefEl_Pk,
                        par :: Parameter.Parameter_MsFEM,
                        problem :: Problem.AbstractPhysicalProblem,
                        k_time :: Int64)
    
    """

    

    """


    
    Phi = [solution.phi_1[:,k_time] solution.phi_2[:,k_time] solution.phi_3[:,k_time]]
    Phi_test = transpose([solution.phi_1[:,k_time] solution.phi_2[:,k_time] solution.phi_3[:,k_time]])
    
    m_loc = Array{Float64,2}(ref_el.n_node, ref_el.n_node)
    
    for k = 1:mesh.n_cell
        m_loc[:,:] = Phi_test * solution.mass * Phi
        for l in 1:9
            M[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += m_loc[l]
        end
    end
    
    return nothing
end


# ---------------------------------------------------------------------------------------------


function assemble_mass_t!(Mt  :: SparseMatrixCSC{Float64,Int64},
                        solution :: Solution_MsFEM,
                        mesh :: Mesh.TriMesh,
                        dof :: FEM.Dof,
                        ref_el :: FEM.RefEl_Pk,
                        par :: Parameter.Parameter_MsFEM,
                        problem :: Problem.AbstractPhysicalProblem,
                        k_time :: Int64)
    
    """

    

    """

    weight_elem = diagm(quad.weight) * Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)
    
    Phi = FEM.eval(ref_el, quad.point)
    Phi_test = FEM.eval(ref_el, quad.point)
    
    m_loc = Array{Float64,2}(ref_el.n_node, ref_el.n_node)
    
    for k = 1:mesh.n_cell
        m_loc[:,:] = Phi_test * diagm(view(weight_elem,:,k)) * Phi'
        for l in 1:9
            M[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += m_loc[l]
        end
    end
    
    return nothing
end


# ---------------------------------------------------------------------------------------------


function assemble_advection!(A  :: SparseMatrixCSC{Float64,Int64},
                             solution :: Solution_MsFEM,
                             mesh :: Mesh.TriMesh,
                             dof :: FEM.Dof,
                             ref_el :: FEM.RefEl_Pk,
                             par :: Parameter.Parameter_MsFEM,
                             problem :: Problem.AbstractPhysicalProblem,
                             k_time :: Int64)

    """

    

    """

    
    return nothing
end


# ---------------------------------------------------------------------------------------------


function assemble_diffusion!(D  :: SparseMatrixCSC{Float64,Int64},
                             solution :: Solution_MsFEM,
                             mesh :: Mesh.TriMesh,
                             dof :: FEM.Dof,
                             ref_el :: FEM.RefEl_Pk,
                             par :: Parameter.Parameter_MsFEM,
                             problem :: Problem.AbstractPhysicalProblem,
                             k_time :: Int64)
    
    """

    

    """


    return nothing
end
