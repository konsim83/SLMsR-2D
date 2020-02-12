function reconstruct_edge_L2_sumEqualOne!(uOpt :: Array{Float64,2},
								 mesh_local  :: Mesh.TriangleMesh.TriMesh, 
								 uGlobal :: Array{Float64,1},
								 uBasis0 :: Array{Float64,2}, 
								 uOrig :: Array{Float64,1}, 
								 ind_seg :: Int,
								 k_edge :: Float64,
								 k_sum :: Float64)

	cell_2d = sort(mesh_local.segment[:,mesh_local.segment_marker.==ind_seg],dims=1)
	if ind_seg==3
		cell_2d[[1;2],end] = cell_2d[[2;1],end]
	end

	# cell_2d = circshift(mesh_local.segment[:,mesh_local.segment_marker.==ind_seg],1)
	ind_edge = (unique(cell_2d))
	n = length(ind_edge)

    basis_lr = zeros(2*n)

    uOrigEdge = uOrig[ind_edge]
    
    if ind_seg==1
	    basis_left = uBasis0[ind_edge,1]
	    basis_right = uBasis0[ind_edge,2]

	    ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = uGlobal[1]
	    a_2 = uGlobal[2]
	elseif ind_seg==2
		basis_left = uBasis0[ind_edge,2]
	    basis_right = uBasis0[ind_edge,3]

	    # ind_con = [1 ; 2 ; 1+n ; 2+n]
	    ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = uGlobal[2] # left
	    a_2 = uGlobal[3] # right
	elseif ind_seg==3
		# Note that here the boundaries of left and right basis are switched
		basis_left = uBasis0[ind_edge,3]
	    basis_right = uBasis0[ind_edge,1]

	    ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = uGlobal[3] # weight of left basis
	    a_2 = uGlobal[1] # weight of right basis
	end

    system_matrix = [ (a_1^2 + k_edge)*speye(n)   (a_1*a_2)*speye(n) ; 
                        (a_1*a_2)*speye(n)   (a_2^2 + k_edge)*speye(n);
                        speye(n)                speye(n)]

    rhs = [ [k_edge*basis_left + a_1*uOrigEdge ;
            k_edge*basis_right + a_2*uOrigEdge ] - system_matrix[1:2*n,ind_con]*val_con ;
            ones(n)]

    basis_lr[ind_con] = val_con
    basis_lr[ind_uncon] = (system_matrix[[ind_uncon ; (2*n+1):3*n], ind_uncon] \ rhs[[ind_uncon ; (2*n+1):3*n]])


    if ind_seg==1
	    uOpt[ind_edge,1] = basis_lr[1:n]
	    uOpt[ind_edge,2] = basis_lr[(n+1):end]
    elseif ind_seg==2
    	uOpt[ind_edge,2] = basis_lr[1:n]
	    uOpt[ind_edge,3] = basis_lr[(n+1):end]
    elseif ind_seg==3
    	# Note that here left and right basis are switched
    	uOpt[ind_edge,3] = basis_lr[1:n]
	    uOpt[ind_edge,1] = basis_lr[(n+1):end]
    end

    return nothing
end
