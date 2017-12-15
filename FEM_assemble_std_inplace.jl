# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------


# ------------------------
# -------   Mass   -------
# ------------------------

function assemble_mass!(M  :: SparseMatrixCSC{Float64,Int64},
                        mesh :: Mesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.RefEl_Pk,
                        quad :: Quad.AbstractQuad,
                        par :: Parameter.Parameter_FEM,
                        problem :: Problem.AbstractProblem)
    
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
# ------------------------
# ------------------------



# -----------------------------
# -------   Advection   -------
# -----------------------------

function assemble_advection!(A  :: SparseMatrixCSC{Float64,Int64},
                             mesh :: Mesh.TriMesh,
                             dof :: FEM.AbstractDof,
                             ref_el :: FEM.RefEl_Pk,
                             quad :: Quad.AbstractQuad,
                             par :: Parameter.Parameter_FEM,
                             problem :: Problem.AbstractProblem,
                             k_time :: Int64)
    
    """

    

    """
    time = k_time * par.dt

    weight_elem = diagm(quad.weight) * Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)
    DF = Mesh.map_ref_point_grad_inv(mesh, quad.point, 1:dof.n_elem);

    x = Mesh.map_ref_point(mesh, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, time, x)
        
    DPhi = FEM.eval_grad(ref_el, quad.point)
    Phi_test = FEM.eval(ref_el, quad.point) 
    
    
    a_loc = zeros(ref_el.n_node,ref_el.n_node)
    
    for k = 1:mesh.n_cell       
        DPhi_mod = modify_ansatzfunction_v_inplace(velocity[:,:,k], DF[:,:,:,k], DPhi)
        a_loc[:,:] = Phi_test * diagm(weight_elem[:,k]) * DPhi_mod'
        for l in 1:9
            A[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += a_loc[l]
        end
    end
    
    return nothing
end


function modify_ansatzfunction_v_inplace(v :: Array{Float64,2}, DF :: Array{Float64,3}, DPhi :: Array{Float64,3})

    DPhi_mod = Array{Float64,2}(size(DPhi,1), size(DPhi,2))

    for i=1:size(DPhi,2)
        DPhi_mod[:,i] = DPhi[:,i,:] * DF[i,:,:] * v[i,:]
    end
    
    return DPhi_mod
end
# -----------------------------
# -----------------------------



# -----------------------------
# -------   Diffusion   -------
# -----------------------------

function assemble_diffusion!(D  :: SparseMatrixCSC{Float64,Int64},
                             m :: Mesh.TriMesh,
                             d :: FEM.AbstractDof,
                             r :: FEM.RefEl_Pk,
                             q :: Quad.AbstractQuad,
                             par :: Parameter.Parameter_FEM,
                             p :: Problem.AbstractProblem,
                             k_time :: Int64)

    n = r.n_node

    time = k_time * par.dt
    
    # fixed quantities for mesh
    weight_elem = diagm(q.weight) * Mesh.map_ref_point_grad_det(m, q.point, 1:d.n_elem)
    DF = Mesh.map_ref_point_grad_inv(m, q.point, 1:d.n_elem)

    x = Mesh.map_ref_point(m, q.point, 1:d.n_elem)
    diffusion = Problem.diffusion(p, time, x)

    DPhi = FEM.eval_grad(r, q.point)
    DPhi_test = FEM.eval_grad(r, q.point)
    
    d_loc = zeros(n,n)

    for ll in 1:m.n_cell
        for kk in 1:length(q.weight)
            tmp = view(DPhi_test,:,kk,:) * view(DF,kk,:,:,ll)
            d_loc += (tmp) * (weight_elem[kk,ll]) *  (tmp * view(diffusion,kk,:,:,ll))'
        end
        
        for l in 1:9
            D[d.ind_test[9*(ll-1)+l],d.ind[9*(ll-1)+l]] -= d_loc[l]
        end
        d_loc .= 0
    end
    
    return nothing
end

# -----------------------------
# -----------------------------
