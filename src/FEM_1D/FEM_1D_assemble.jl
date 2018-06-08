# !!! Mass 1D !!!
function assemble_mass_1d(mesh :: Mesh_1D,
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
    
    mat_local = assemble_elem_r(mesh, quad, ref_el, problem)

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
    
    
    Phi = shapefun(ref_el, quad.point) * diagm(quad.weight)
    Phi_x = shapefun_der(ref_el, quad.point)

    

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
            DtVel = quad.weight[j] * velocity[k][j]
            for i2=1:ref_el.n_basis
                for i1=1:ref_el.n_basis
                    r[i1,i2,k] -= (Phi[i1,j]  * dot(DtVel,(p_end[:,i] - p_start[:,i])/norm(p_end[:,i] - p_start[:,i])) 
                                    + dot(DtVel,DPhi[i1,j])  * Phi[i2,j] )
                    # a[i1,i2,k] += Phi_test[i1,j]  * dot(DtVel,DPhi[i2,j])
                end
            end            
        end
    end

    return r
end # end function


#=
# ----------------------------------------------------------



function assemble_diffusion_1d(cells :: Array{Int,2},
                            points :: Array{Int,2},
                            ref_el :: RefEl_Lagrange_1,
                            quad :: Quad1D.JacobiGauss,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    ind = vec(transpose([cells[1,:] cells[1,:] cells[2,:] cells[2,:]]))
    ind_test = vec([cells ; cells])

    p_start = points[:,2,]
    p_end = points[:,1,]
    
    mat_local = assemble_elem_d(p_start, p_end, ref_el , quad, problem, t)

    Mat_global = sparse(ind_test, ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end # end function

function assemble_elem_d(cells :: Array{Int,2},
                            p_start :: Array{Int,2},
                            p_end :: Array{Int,2},
                            ref_el :: AbstractRefEl_Pk,
                            quad :: Quad1D.JacobiGauss,
                            problem :: Problem.AbstractProblem,
                            t :: Float64)

    n = length(ref_el.node)
    d = zeros(n, n, size(cell,2))
    
    # x must be an array of 2-by-n_q arrays
    x = [hcat([(p_end[:,i] - p_start[:,i])*q + p_start[:,i] for q in quad.point]...) for i in 1:size(cell,2)]
    diffusion = Problem.diffusion(problem, t, x)
    
    weight_elem = reshape( sqrt.(sum((p_end-p_start).^2,1)) , 1, 1, size(cells,2))


    Phi = eval_shapefun(ref_el, quad.point) * diagm(quad.weight)
    Phi_x = eval_shapefun_der(ref_el, quad.point)

    for k = 1:mesh.n_cell
        for j=1:quad.n_point
                DtVel =  quad.weight[j] * diffusion[k][j] / weight_elem[k]
                for i2=1:ref_el.n_basis
                    for i1=1:ref_el.n_basis
                        d[i1,i2,k] += dot( DPhi_test[i1,j] ,  DtDiffD*DPhi[i2,j])
                        # a[i1,i2,k] += Phi_test[i1,j]  * dot(DtVel,DPhi[i2,j])
                    end
                end
            end
        end
    end

    return d
end # end function
=#