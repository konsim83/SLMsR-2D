# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------


# ------------------------
# -------   Mass   -------
# ------------------------

function assemble_mass!(M  :: SparseMatrixCSC{Float64,Int64},
                        mesh :: Mesh.TriangleMesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.RefEl_Pk{N},
                        quad :: Quad.AbstractQuad,
                        problem :: Problem.AbstractProblem) where N

    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    Phi = FEM.shapeFun(ref_el, quad.point)
    Phi_test = FEM.shapeFun(ref_el, quad.point) * diagm(quad.weight)
    
    m_loc = Array{Float64,2}(ref_el.n_basis, ref_el.n_basis)

    Nsquare = ref_el.n_basis^2
    for k = 1:mesh.n_cell
        m_loc[:,:] = Phi_test * weight_elem[k] * Phi'
        for l in 1:Nsquare
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
                             mesh :: Mesh.TriangleMesh.TriMesh,
                             dof :: FEM.AbstractDof,
                             ref_el :: FEM.RefEl_Pk{N},
                             quad :: Quad.AbstractQuad,
                             par :: Parameter.AbstractParameter,
                             problem :: Problem.AbstractProblem,
                             k_time :: Int64) where N
    
    time_idx = k_time * par.dt 
        
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, time_idx, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem)
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    DPhi = FEM.shapeFun_grad(ref_el, quad.point)
    Phi_test = FEM.shapeFun(ref_el, quad.point)

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
            DtVel = weight_elem[k] * quad.weight[j] * DF[k] * velocity[k][j]
            l = 0
            for i2=1:ref_el.n_basis
                for i1=1:ref_el.n_basis
                    l += 1
                    A[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] += Phi_test[i1,j]  * dot(DtVel,DPhi[i2,j])
                end
            end
        end
    end

    
    return nothing
end
# -----------------------------
# -----------------------------



# -----------------------------
# -------   Diffusion   -------
# -----------------------------
function assemble_diffusion!(D  :: SparseMatrixCSC{Float64,Int64},
                             mesh :: Mesh.TriangleMesh.TriMesh,
                             dof :: FEM.AbstractDof,
                             ref_el :: FEM.RefEl_Pk{N},
                             quad :: Quad.AbstractQuad,
                             par :: Parameter.AbstractParameter,
                             problem :: Problem.AbstractProblem,
                             k_time :: Int64) where N
    
    time_idx = k_time * par.dt 
        
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    diffusion = Problem.diffusion(problem, time_idx, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem);
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    DPhi = FEM.shapeFun_grad(ref_el, quad.point)
    DPhi_test = FEM.shapeFun_grad(ref_el, quad.point)

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
            DtDiffD = weight_elem[k] * quad.weight[j] * DF[k]*diffusion[k][j]*DF[k]'
            l = 0
            for i2=1:ref_el.n_basis
                for i1=1:ref_el.n_basis
                    l += 1
                    D[dof.ind_test[9*(k-1)+l],dof.ind[9*(k-1)+l]] -= dot( DPhi_test[i1,j] , DtDiffD*DPhi[i2,j])
                end
            end
        end
    end

    
    return nothing
end
# -----------------------------
# -----------------------------
