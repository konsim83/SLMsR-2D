function reconstruct_edge!(uOpt, m_f, U, u0, u, ind_seg, k_edge)

	cell_2d = circshift(m_f.segment[:,m_f.segment_marker.==ind_seg],1)
	ind_edge = sort(unique(cell_2d))
	n = length(ind_edge)

    basis_lr = zeros(2*n)

    u_before = u[ind_edge]
    
    if ind_seg==1
	    basis_left = u0[ind_edge,1]
	    basis_right = u0[ind_edge,2]

	    ind_con = [1 ; 2 ; 1+n ; 2+n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[1]
	    a_2 = U[2]
	elseif ind_seg==2
		basis_left = u0[ind_edge,2]
	    basis_right = u0[ind_edge,3]

	    ind_con = [1 ; 2 ; 1+n ; 2+n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[2]
	    a_2 = U[3]
	elseif ind_seg==3
		basis_left = u0[ind_edge,1]
	    basis_right = u0[ind_edge,3]

	    ind_con = [1 ; 2 ; 1+n ; 2+n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[1]
	    a_2 = U[3]
	end

    system_matrix = [ (a_1^2 + k_edge)*speye(n)   (a_1*a_2)*speye(n) ; 
                        (a_1*a_2)*speye(n)   (a_2^2 + k_edge)*speye(n) ]

    rhs = ( [k_edge*basis_left + a_1*u_before ;
            k_edge*basis_right + a_2*u_before]
            - system_matrix[:,ind_con]*val_con )
    
    basis_lr[ind_con] = val_con
    basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]


    if ind_seg==1
	    uOpt[ind_edge,1] = basis_lr[1:n]
	    uOpt[ind_edge,2] = basis_lr[(n+1):end]
    elseif ind_seg==2
    	uOpt[ind_edge,2] = basis_lr[1:n]
	    uOpt[ind_edge,3] = basis_lr[(n+1):end]
    elseif ind_seg==3
    	# Note that here left and right basis are switched
    	uOpt[ind_edge,3] = basis_lr[(n+1):end]
	    uOpt[ind_edge,1] = basis_lr[1:n]
    end

    return nothing
end