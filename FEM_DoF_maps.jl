# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

function map_ref_point(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :

    P = [dof.T_ref2cell[ind] * [x  ; ones(1,size(x,2))] for ind in ind_c]
    
    return P
end


function map_ref_point(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    return map_ref_point(dof, x, collect(ind_c))
end


function map_ref_point(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Int64)

    return map_ref_point(dof, x, [ind_c])[1]
end


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


function map_ref_point_grad(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :
    
    n_points = size(x,1)
    n_cells = length(ind_c)
    
    P_grad = [dof.T_ref2cell[ind][:,1:2] for ind in ind_c]
    
    return P_grad
end


function map_ref_point_grad(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    return map_ref_point_grad(dof, x, collect(ind_c))
end


function map_ref_point_grad(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Int64)

    return map_ref_point_grad(dof, x, [ind_c])[1]
end


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


function map_ref_point_grad_inv(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :
    
    n_points = size(x,1)
    n_cells = length(ind_c)
    
    P_grad_inv = [dof.T_cell2ref[ind][:,1:2] for ind in ind_c]

    return P_grad_inv
end


function map_ref_point_grad_inv(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})

    return map_ref_point_grad_inv(dof, x, collect(ind_c))
end


function map_ref_point_grad_inv(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Int64)

    return map_ref_point_grad_inv(dof, x, [ind_c])[1]
end


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


function map_ref_point_grad_det(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Array{Int64,1})
    
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :
    
    P_det = [abs(det(dof.T_ref2cell[ind][1:2,1:2])) for ind in ind_c]

    return P_det
end


function map_ref_point_grad_det(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: UnitRange{Int64})
    
    return map_ref_point_grad_det(dof, x, collect(ind_c))
end


function map_ref_point_grad_det(dof :: AbstractDof, x :: Array{Float64,2}, ind_c :: Int64)
    
    return map_ref_point_grad_det(dof, x, [ind_c])[1]
end


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------