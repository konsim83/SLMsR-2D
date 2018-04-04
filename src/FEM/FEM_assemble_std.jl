# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------

# ------------------------
# -------   Mass   -------
# ------------------------
function assemble_mass(mesh :: Mesh.TriangleMesh.TriMesh,
                       dof :: FEM.AbstractDof,
                       ref_el :: FEM.RefEl_Pk{N},
                       quad :: Quad.AbstractQuad,
                       problem :: Problem.AbstractProblem) where N
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_m(mesh, ref_el, dof, quad)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end


function assemble_elem_m(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk{N},
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad) where N
    
    n = ref_el.n_basis
    mass = Array{Float64,3}(n, n, mesh.n_cell)
    
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    Phi = FEM.shapeFun(ref_el, quad.point)
    Phi_test = FEM.shapeFun(ref_el, quad.point) * diagm(quad.weight)

    for k = 1:mesh.n_cell
        mass[:,:,k] = Phi_test * weight_elem[k] * Phi'
    end

    return mass
end # end function



# -----------------------------
# -------   Advection   -------
# -----------------------------
function assemble_advection(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            problem :: Problem.AbstractProblem,
                            time_idx :: Float64)
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_a(mesh, ref_el, dof, quad, problem, time_idx)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_a(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time_idx :: Float64)
    
    n = ref_el.n_basis
    a = zeros(n, n, mesh.n_cell)
    
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, time_idx, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem)
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    if problem.conservative
        DPhi_test = FEM.shapeFun_grad(ref_el, quad.point)
        Phi = FEM.shapeFun(ref_el, quad.point)

        for k = 1:mesh.n_cell
            for j=1:quad.n_point
                DtVel = weight_elem[k] * quad.weight[j] * DF[k] * velocity[k][j]
                for i2=1:ref_el.n_basis
                    for i1=1:ref_el.n_basis
                        a[i1,i2,k] -= dot(DtVel,DPhi_test[i1,j])  * Phi[i2,j]
                    end
                end
            end
        end
    else
        DPhi = FEM.shapeFun_grad(ref_el, quad.point)
        Phi_test = FEM.shapeFun(ref_el, quad.point)

        for k = 1:mesh.n_cell
            for j=1:quad.n_point
                DtVel = weight_elem[k] * quad.weight[j] * DF[k] * velocity[k][j]
                for i2=1:ref_el.n_basis
                    for i1=1:ref_el.n_basis
                        a[i1,i2,k] += Phi_test[i1,j]  * dot(DtVel,DPhi[i2,j])
                    end
                end
            end
        end
    end

    return a
end # end function



# -----------------------------
# -------   Diffusion   -------
# -----------------------------
function assemble_diffusion(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            problem :: Problem.AbstractProblem,
                            time_idx :: Float64)
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_d(mesh, ref_el, dof, quad, problem, time_idx)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_d(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time_idx :: Float64)
    
    n = ref_el.n_basis
    d_loc = zeros(n, n, mesh.n_cell)
    
    
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    diffusion = Problem.diffusion(problem, time_idx, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem);
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    DPhi = FEM.shapeFun_grad(ref_el, quad.point)
    DPhi_test = FEM.shapeFun_grad(ref_el, quad.point)

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
            DtDiffD = weight_elem[k] * quad.weight[j] * DF[k]*diffusion[k][j]*DF[k]'
            for i2=1:ref_el.n_basis
                for i1=1:ref_el.n_basis
                    # d_loc[i1,i2,k] = weight_elem[k] * quad.weight[j] * (DPhi_test[i1,j]'*DF[k]) * (diffusion[k][j]*(DF[k]'*DPhi[i2,j]))
                    d_loc[i1,i2,k] += dot( DPhi_test[i1,j] ,  DtDiffD*DPhi[i2,j])
                end
            end
        end
    end

    return -d_loc
end # end function




# -----------------------------
# -------   Advection   -------
# -----------------------------
function assemble_reaction(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            problem :: Problem.AbstractProblem,
                            time_idx :: Float64)
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_r(mesh, ref_el, dof, quad, problem, time_idx)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_r(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time_idx :: Float64)
    
    n = ref_el.n_basis
    r = zeros(n, n, mesh.n_cell)
    
    if problem.conservative
        x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
        velocity = Problem.velocity(problem, time_idx, x)

        DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem)
        weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

        DPhi = FEM.shapeFun_grad(ref_el, quad.point)
        Phi = FEM.shapeFun(ref_el, quad.point)
    
        for k = 1:mesh.n_cell
            for j=1:quad.n_point
                DtVel = weight_elem[k] * quad.weight[j] * DF[k] * velocity[k][j]
                for i2=1:ref_el.n_basis
                    for i1=1:ref_el.n_basis
                        r[i1,i2,k] -= Phi[i1,j]  * dot(DtVel,DPhi[i2,j]) + dot(DtVel,DPhi[i1,j])  * Phi[i2,j]
                    end
                end
            end
        end
    end

    return r
end # end function

# -----------------------------
# -----------------------------