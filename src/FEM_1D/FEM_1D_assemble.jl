function assemble_mass(mesh :: Mesh_1D,
                        quad :: Quad.Quad_line,
                        ref_el :: RefEl_Lagrange_1,
                        dof :: Dof_1D,
                        problem :: Problem.AbstractProblem)

    ind_trial = dof.ind_trial
    ind_test = dof.ind_test
    
    weight_elem = map_ref_point_der(mesh, quad.point, collect(1:mesh.n_cell))
    mat_local = shapeFun(ref_el, quad.point) * diagm(quad.weight) * transpose(shapeFun(ref_el, quad.point))
    
    mat_local = hcat([mat_local*w for w in weight_elem]...)

    Mat_global = sparse(ind_test, ind_trial, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end # end function



# ----------------------------------------------------------



function assemble_reaction(mesh :: Mesh_1D,
                            quad :: Quad.Quad_line,
                            ref_el :: RefEl_Lagrange_1,
                            dof :: Dof_1D,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    ind_trial = dof.ind_trial
    ind_test = dof.ind_test
    
    mat_local = assemble_elem_r(mesh, quad, ref_el, problem, t)

    Mat_global = sparse(ind_test, ind_trial, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end # end function

function assemble_elem_r(mesh :: Mesh_1D,
                            quad :: Quad.Quad_line,
                            ref_el :: RefEl_Lagrange_1,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    n = length(ref_el.node)
    r = zeros(n, n, mesh.n_cell)
    
    # x must be an array of 2-by-n_q arrays
    x = map_ref_point(mesh, quad.point, collect(1:mesh.n_cell))
    velocity = Problem.velocity(problem, t, x)
    normal_elem = get_cell_normal(mesh, collect(1:mesh.n_cell))    
    
    Phi = shapeFun(ref_el, quad.point) * diagm(quad.weight)
    Phi_x = shapeFunDer(ref_el, quad.point)    

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
            Vel = dot(normal_elem[k],velocity[k][j])
            r[:,:,k] -= Phi[:,j] * Vel * transpose(Phi_x[:,j]) + Phi_x[:,j] * Vel * transpose(Phi[:,j])
        end
    end

    return r
end # end function



# ----------------------------------------------------------



function assemble_diffusion(mesh :: Mesh_1D,
                            quad :: Quad.Quad_line,
                            ref_el :: RefEl_Lagrange_1,
                            dof :: Dof_1D,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    ind_trial = dof.ind_trial
    ind_test = dof.ind_test
    
    mat_local = assemble_elem_d(mesh, quad, ref_el, problem, t)

    Mat_global = sparse(ind_test, ind_trial, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end # end function


function assemble_elem_d(mesh :: Mesh_1D,
                            quad :: Quad.Quad_line,
                            ref_el :: RefEl_Lagrange_1,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    n = length(ref_el.node)
    d = zeros(n, n, mesh.n_cell)
    
    # x must be an array of 2-by-n_q arrays
    x = map_ref_point(mesh, quad.point, collect(1:mesh.n_cell))
    diffusion = Problem.diffusion(problem, t, x)
    normal_elem = get_cell_normal(mesh, collect(1:mesh.n_cell))
    weight_elem = map_ref_point_der(mesh, quad.point, collect(1:mesh.n_cell))
    
    Phi_x = shapeFunDer(ref_el, quad.point) * diagm(quad.weight)

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
            Dif = transpose(normal_elem[k]) * diffusion[k][j] * normal_elem[k]
            d[:,:,k] += Phi_x[:,j] * Dif * transpose(Phi_x[:,j]) / weight_elem[k]
        end
    end

    return -d
end # end function