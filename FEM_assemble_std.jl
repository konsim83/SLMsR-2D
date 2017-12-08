# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------


# ------------------------
# -------   Mass   -------
# ------------------------
function assemble_mass(mesh :: Mesh.TriMesh,
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
    mat_local = assemble_elem_m(mesh, ref_el, dof, quad, problem, time)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end


function assemble_elem_m(mesh :: Mesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time :: Float64)
    
    n = ref_el.n_node
    m = Array{Float64,3}(n, n, mesh.n_cell)
    
    # x = Mesh.map_ref_point(mesh, quad.point, 1:dof.n_elem)
    
    weight_elem = Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)

    Phi = FEM.eval(ref_el, quad.point)
    Phi_test = FEM.eval(ref_el, quad.point)

    for k = 1:mesh.n_cell    
        m[:,:,k] = Phi_test * diagm(quad.weight) * diagm(weight_elem[:,k]) * Phi'
    end

    return m
end # end function


# ------------------------
# ------------------------



# -----------------------------
# -------   Advection   -------
# -----------------------------
function assemble_advection(mesh :: Mesh.TriMesh,
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


function assemble_elem_a(mesh :: Mesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time :: Float64)
    
    n = ref_el.n_node
    a = zeros(n, n, mesh.n_cell)    
    
    x = Mesh.map_ref_point(mesh, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, time, x)

    DF = Mesh.map_ref_point_grad_inv(mesh, quad.point, 1:dof.n_elem);
    weight_elem = Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)

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
function assemble_diffusion(mesh :: Mesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.AbstractProblem,
                            k_time :: Int64)

    """

    

    """

    # Assembly pattern
    #=
    i = get_dof_elem(dof, mesh, 1:dof.n_elem)
    ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
    ind_test = vec(transpose(repmat(i, 1, size(i,2))))
    =#

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_d(mesh, ref_el, dof, quad, problem, time)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_d(mesh :: Mesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time :: Float64)
    
    n = ref_el.n_node
    a = zeros(n, n, mesh.n_cell)
    
    
    x = Mesh.map_ref_point(mesh, quad.point, 1:dof.n_elem)
    diffusion = Problem.diffusion(problem, time, x)

    DF = Mesh.map_ref_point_grad_inv(mesh, quad.point, 1:dof.n_elem);
    weight_elem = Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)

    Phi = FEM.eval_grad(ref_el, quad.point)
    Phi_test = FEM.eval_grad(ref_el, quad.point)
    
    for k = 1:mesh.n_cell    
        a[:,:,k] = build_elem_matrix_d(diffusion[:,:,:,k], DF[:,:,:,k], Phi_test, Phi, quad.weight, weight_elem[:,k])
    end

    return a
end # end function

function build_elem_matrix_d(diff :: Array{Float64,3}, DF :: Array{Float64,3},
                                 Phi_test :: Array{Float64,3}, Phi :: Array{Float64,3},
                                 q_weight :: Array{Float64,1}, weight_elem :: Array{Float64,1})

    Phi_mod = zeros(size(Phi,1), size(Phi,2))

    for i=1:size(Phi,1)
        Phi_mod += (Phi_test[i,:,:] * DF[i,:,:]) * (q_weight[i] * weight_elem[i]) *  (Phi[i,:,:] * DF[i,:,:] * diff[i,:,:])'
    end
    
    return -Phi_mod
end

# -----------------------------
# -----------------------------
