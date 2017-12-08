function map_ref_point(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    y = [x  ones(size(x,1))]
    n_points = size(x,1)

    T = reshape(permutedims(m.T_ref2cell[:,:,ind_c],[2,1,3]),3,2*length(ind_c))
    P = reshape(y*T, n_points, 2, length(ind_c))
    
    return P
end


function map_ref_point(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    ind_c = vec(collect(ind_c))
    
    return map_ref_point(m, x, ind_c)
end


# --------------------------------------------------------------------------


function map_ref_point_grad(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)
    
    P_grad = Array{Float64,4}(n_points, 2, 2, n_cells)

    for j=1:n_points
        P_grad[j,:,:,:] = m.T_ref2cell[:,1:2,ind_c]
    end
    
    return P_grad
end


function map_ref_point_grad(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    ind_c = vec(collect(ind_c))
    
    return map_ref_point_grad_inv(m, x, ind_c)
end


# --------------------------------------------------------------------------


function map_ref_point_grad_inv(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    ind_c = vec(collect(ind_c))
    
    n_points = size(x,1)
    n_cells = length(ind_c)
    
    P_grad = Array{Float64,4}(n_points, 2, 2, n_cells)

    for j=1:n_points
        P_grad[j,:,:,:] = m.T_cell2ref[:,1:2,ind_c]
    end
    
    return P_grad
end


function map_ref_point_grad_inv(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    ind_c = vec(collect(ind_c))
    
    return map_ref_point_grad_inv(m, x, ind_c)
end


# --------------------------------------------------------------------------


function map_ref_point_grad_det(m :: TriMesh, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    n_points = size(x,1)
    n_cells = length(ind_c)

    P_det = Array{Float64,2}(n_points, n_cells)
    
    for i=1:n_cells
        for j=1:n_points
            P_det[j,i] = abs(m.T_ref2cell[1,1,ind_c[i]]*m.T_ref2cell[2,2,ind_c[i]]-m.T_ref2cell[2,1,ind_c[i]]*m.T_ref2cell[1,2,ind_c[i]])
        end
    end
    
    return P_det
end


function map_ref_point_grad_det(m :: TriMesh, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})
    ind_c = vec(collect(ind_c))
    
    return map_ref_point_grad_det(m, x, ind_c)
end


# --------------------------------------------------------------------------
