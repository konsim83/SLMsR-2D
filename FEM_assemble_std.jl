# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------


# ------------------------
# -------   Mass   -------
# ------------------------
function assemble_mass(mesh :: Mesh.TriangleMesh.TriMesh,
                       dof :: FEM.AbstractDof,
                       ref_el :: FEM.RefEl_Pk,
                       quad :: Quad.AbstractQuad,
                       par :: Parameter.Parameter_FEM,
                       problem :: Problem.AbstractProblem)
    
    """

    

    """
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_m(mesh, ref_el, dof, quad)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end


function assemble_elem_m(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad)
    
    n = ref_el.n_node
    m = Array{Float64,3}(n, n, mesh.n_cell)
    
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    Phi = FEM.eval(ref_el, quad.point)
    Phi_test = FEM.eval(ref_el, quad.point) * diagm(quad.weight)

    for k = 1:mesh.n_cell    
        m[:,:,k] = Phi_test  * diagm(weight_elem[:,k]) * Phi'
    end

    return m
end # end function


# ------------------------
# ------------------------



# -----------------------------
# -------   Advection   -------
# -----------------------------
function assemble_advection(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.AbstractProblem,
                            k_time :: Int64)
    
    """

    

    """

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_a(mesh, ref_el, dof, quad, problem, time)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_a(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time :: Float64)
    
    n = ref_el.n_node
    a = zeros(n, n, mesh.n_cell)    
    
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, time, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem);
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    Phi = FEM.eval_grad(ref_el, quad.point)
    Phi_test = FEM.eval(ref_el, quad.point)
    
    for k = 1:mesh.n_cell    
        a[:,:,k] = Phi_test * diagm(quad.weight) * diagm(weight_elem[:,k]) * modify_ansatzfunction_v(velocity[:,:,k], DF[:,:,:,k], Phi)'
    end

    return a
end # end function

function modify_ansatzfunction_v(v :: Array{Float64,2}, DF :: Array{Float64,3}, DPhi :: Array{Float64,3})

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
function assemble_diffusion(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.AbstractProblem,
                            k_time :: Int64)

    """

    

    """

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_d(mesh, ref_el, dof, quad, problem, time)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_d(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time :: Float64)
    
    n = ref_el.n_node
    d_loc = zeros(n, n, mesh.n_cell)
    
    
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    diffusion = Problem.diffusion(problem, time, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem);
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    DPhi = FEM.eval_grad(ref_el, quad.point)
    DPhi_test = FEM.eval_grad(ref_el, quad.point)
    
    for l = 1:mesh.n_cell
        for k = 1:length(quad.weight)
            d_loc[:,:,l] += (view(DPhi_test,:,k,:) * view(DF,k,:,:,l)) * (quad.weight[k] * weight_elem[k,l]) *  (view(DPhi,:,k,:) * view(DF,k,:,:,l) * view(diffusion,k,:,:,l))'
        end
        #d[:,:,k] = build_elem_matrix_d(diffusion[:,:,:,k], DF[:,:,:,k], DPhi_test, DPhi, quad.weight, weight_elem[:,k])
    end

    return -d_loc
end # end function

function build_elem_matrix_d(diff :: Array{Float64,3}, DF :: Array{Float64,3},
                                 DPhi_test :: Array{Float64,3}, DPhi :: Array{Float64,3},
                                 q_weight :: Array{Float64,1}, weight_elem :: Array{Float64,1})

    DPhi_mod = zeros(size(DPhi,1), size(DPhi,2))
    
    for k = 1:length(q_weight)        
        DPhi_mod += (DPhi_test[:,k,:] * DF[k,:,:]) * (q_weight[k] * weight_elem[k]) *  (DPhi[:,k,:] * DF[k,:,:] * diff[k,:,:])'
    end
    
    return -DPhi_mod
end

# -----------------------------
# -----------------------------
