# --------------------   Reaction term in Lagrangian coordinnates   --------------------
# --------------------   u_t = Du - <c_x u, v> = Du + <cu_x, v> + <cu, v_x> = (D+R)u
function assemble_mass_1(cells :: Array{Int,2},
                            ref_el :: RefEl_Lagrange_1,
                            quad :: Quad1D.JacobiGauss,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    ind = get_dofs(dof, collect(1:dof.n_elem))
    ind_test = vec(transpose(repmat(ind, 1, size(ind,2))))
    ind = vec(sortrows(transpose(repmat(ind, 1, size(ind,2)))))

    mat_local = assemble_elem_m(mesh, ref_el , quad, problem, t)

    Mat_global = sparse(ind_test, ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end # end function


function assemble_elem_r(mesh :: Grid1D.Mesh,
                            ref_el :: AbstractRefEl_Pk,
                            quad :: Quad1D.JacobiGauss,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    n = length(ref_el.node)
    r = zeros(n, n, mesh.n_cell)
    
    x = Grid1D.map_ref_point(mesh, quad.point, collect(1:mesh.n_cell))
    y = x #get_domain_points(x, [mesh.a; mesh.b])

    velocity = Problem.velocity(problem, t, y)
    weight_elem = Grid1D.map_ref_point_der(mesh, 0.5, collect(1:mesh.n_cell))

    Phi = eval_shapefun(ref_el, quad.point) * diagm(quad.weight)
    Phi_x = eval_shapefun_der(ref_el, quad.point)

    for k = 1:mesh.n_cell
        r[:,:,k] = build_local_reaction_matrix(velocity[k,:], 
                                                weight_elem[k],
                                                Phi, Phi_x)
    end

    return r
end # end function


function build_local_reaction_matrix(velocity :: Array{Float64,1}, 
                                        weight_elem :: Float64,
                                        Phi :: Array{Float64,2},
                                        Phi_x :: Array{Float64,2})
        
    V = repmat(velocity', size(Phi,1), 1)

    # u_t + <c_x u, v> = u_t - <cu_x, v> - <cu, v_x> = u_t + Ru = Du
    r = - ( (Phi) * transpose(Phi_x.*V) 
            + (Phi_x) * transpose(Phi.*V) )

    return r
end # end function
