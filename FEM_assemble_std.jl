# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------


# ------------------------
# -------   Mass   -------
# ------------------------

function assemble_mass(solution :: Solution_FEM,
                       mesh :: Mesh.TriMesh,
                       dof :: FEM.Dof,
                       ref_el :: FEM.RefEl_Pk,
                       quad :: Quadrature.Quad_simplex,
                       par :: Parameter.Parameter_FEM,
                       problem :: Problem.Problem_top,
                       k_time :: Int64)
    
    """

    

    """

    # Assembly pattern is []
    i = get_dof_elem(dof, mesh, 1:dof.n_elem)
    ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
    ind_test = vec(transpose(repmat(i, 1, size(i,2))))

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_m(mesh, ref_el , quad, problem, time)

    Mat_global = sparse(ind_test, ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end


function assemble_elem_m(mesh :: Mesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         quad :: Quadrature.Quad_simplex,
                         problem :: Problem.Problem_top,
                         time :: Float64)
    
    n = length(ref_el.n_node)
    m = zeros(n, n, mesh.n_elem)
    
    # x = Mesh.map_ref_point(mesh, quad.point, 1:dof.n_elem)
    
    weight_elem = Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)

    Phi = FEM.eval(r, q.point)
    Phi_test = FEM.eval(r, q.point)

    for k = 1:mesh.n_elem    
        m[:,:,k] = Phi_test' * diagm(quad.weight) * diagm(weight_elem[:,k]) * Phi
    end

    return m
end # end function


# ------------------------
# ------------------------



# -----------------------------
# -------   Advection   -------
# -----------------------------

function assemble_advection(solution :: Solution_FEM,
                            mesh :: Mesh.TriMesh,
                            dof :: FEM.Dof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quadrature.Quad_simplex,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.Problem_top,
                            k_time :: Int64)
    
    """

    

    """

    # Assembly pattern is []
    i = d.get_dof_elem(1:d.n_elem)
    ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
    ind_test = vec(transpose(repmat(i, 1, size(i,2))))

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_a(mesh, ref_el , quad, problem, time)

    Mat_global = sparse(ind_test, ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_a(mesh :: Mesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         quad :: Quadrature.Quad_simplex,
                         problem :: Problem.Problem_top,
                         time :: Float64)
    
    n = length(ref_el.n_node)
    a = zeros(n, n, mesh.n_elem)
    
    
    x = Mesh.map_ref_point(mesh, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, t, x)
    
    weight_elem = Mesh.map_ref_point_grad_det(mesh, quad.point, 1:dof.n_elem)

    Phi = FEM.eval_der(r, q.point)
    Phi_test = FEM.eval(r, q.point)

    for k = 1:mesh.n_elem    
        a[:,:,k] = Phi_test' * diagm(quad.weight) * diagm(weight_elem[:,k]) * Phi
    end

    return a
end # end function

# -----------------------------
# -----------------------------




# -----------------------------
# -------   Diffusion   -------
# -----------------------------

function assemble_diffusion(solution :: Solution_FEM,
                            mesh :: Mesh.TriMesh,
                            dof :: FEM.Dof,
                            ref_el :: FEM.RefEl_Pk,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.Problem_top,
                            k_time :: Int64)

    """

    

    """

end

# -----------------------------
# -----------------------------
