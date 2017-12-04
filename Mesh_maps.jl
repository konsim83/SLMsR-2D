function map_ref_point(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    y = transpose( [x  ones(size(x,1))] )
    n_points = size(x,1)
    n_cells = length(ind_c)

    s = "sparse($y[:,:]),"^n_cells
    s1 = s[1:end-1]
    s2 = string("blkdiag(", s1, ")")
    e = parse(s2)
    D = eval(e)

    s = ""
    for i=1:n_cells
        T = m.T_ref2cell[:,:,ind_c[i]]
        s = string(s, "$T,")
    end
    s1 = s[1:end-1]
    s2 = string("hcat(", s1, ")")
    e = parse(s2)
    A = eval(e)

    P0 = full(A*D)
    P0 = reshape(P0, 2, n_points, n_cells)
    P1 = reshape(P0[1,:], n_points, 1, n_cells)
    P2 = reshape(P0[2,:], n_points, 1, n_cells)
    P = hcat(P1, P2)
    
    return P
end


function map_ref_point(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    ind_c = vec(collect(ind_c))

    y = transpose( [x  ones(size(x,1))] )
    n_points = size(x,1)
    n_cells = length(ind_c)

    s = "sparse($y[:,:]),"^n_cells
    s1 = s[1:end-1]
    s2 = string("blkdiag(", s1, ")")
    e = parse(s2)
    D = eval(e)

    s = ""
    for i=1:n_cells
        T = m.T_ref2cell[:,:,ind_c[i]]
        s = string(s, "$T,")
    end
    s1 = s[1:end-1]
    s2 = string("hcat(", s1, ")")
    e = parse(s2)
    A = eval(e)

    P0 = full(A*D)
    P0 = reshape(P0, 2, n_points, n_cells)
    P1 = reshape(P0[1,:], n_points, 1, n_cells)
    P2 = reshape(P0[2,:], n_points, 1, n_cells)
    P = hcat(P1, P2)
    
    return P
end


# --------------------------------------------------------------------------


function map_ref_point_grad(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_grad = Array{Float64,4}(n_points, 2, 2, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_grad[j,:,:,i] = m.T_ref2cell[:,1:2,ind_c[i]]
        end
    end
    
    return P_grad
end


function map_ref_point_grad(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_grad = Array{Float64,4}(n_points, 2, 2, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_grad[j,:,:,i] = m.T_ref2cell[:,1:2,ind_c[i]]
        end
    end
    
    return P_grad
end


# --------------------------------------------------------------------------


function map_ref_point_grad_inv(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_grad = Array{Float64,4}(n_points, 2, 2, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_grad[j,:,:,i] = m.T_cell2ref[:,1:2,ind_c[i]]
        end
    end
    
    return P_grad
end


function map_ref_point_grad_inv(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_grad = Array{Float64,4}(n_points, 2, 2, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_grad[j,:,:,i] = m.T_cell2ref[:,1:2,ind_c[i]]
        end
    end
    
    return P_grad
end


# --------------------------------------------------------------------------


function map_ref_point_grad_det(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_det = Array{Float64,2}(n_points, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_det[j,i] = abs(det(m.T_ref2cell[:,1:2,ind_c[i]]))
        end
    end
    
    return P_det
end


function map_ref_point_grad_det(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_det = Array{Float64,2}(n_points, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_det[j,i] = abs(det(m.T_ref2cell[:,1:2,ind_c[i]]))
        end
    end
    
    return P_det
end


# --------------------------------------------------------------------------
