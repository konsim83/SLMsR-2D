struct RefEl_Pk{N} <: AbstractRefEl
    info_str :: String

    n_order :: Int64
    
    node :: Array{Float64, 2}
    n_node :: Int64

    coeff :: Array{Float64, 2}
    
end # end type


function RefEl_Pk(n_order :: Int)

    info_str = "Triangular Lagrange element of type P$n_order."

    if n_order==1
        # -------   P1   -------
        n_order = 1
        
        node = [0 1 0 ; 0 0 1]
        n_node = 3

        # columns are coefficients of basis functions
        # phi(x) = ax + by + c
        coeff = [ node ; 
                    ones(1,n_node) ] \ eye(n_node)

    elseif n_order==2
            # -------   P2   -------
            n_order = 2
            
            node = [0 1 0 0.5 0.5 0 ; 0 0 1 0 0.5 0.5]
            n_node = 6

            # columns are coefficients of basis functions
            # phi(x) = ax^2 + by^2 + cxy + dx + ey + f
            coeff = [ node.^2 ; 
                        node[1,:]'.*node[2,:]' ; 
                        node ; 
                        ones(1, n_node) ] \ eye(n_node)

    elseif n_order==3
            # -------   P3   -------
            n_order = 3
            
            node = [0 1 0 1/3 2/3 2/3 1/3 0 0 1/2 ; 0 0 1 0 0 1/3 2/3 2/6 1/3 1/2]
            n_node = 10

            # columns are coefficients of basis functions
            # phi(x) = ax^3 + by^3 + cx^2y + dxy^2 + ex^2 + fy^2 + gxy + hx + iy + j
            coeff = [ node.^3 ; node[1,:]'.^2.*node[2,:]' ; node[1,:]'.*node[2,:]'.^2 ; 
                        node[1,:].^2 ; node[2,:].^2 ; node[1,:]'.*node[2,:]' ; 
                        node ; 
                        ones(1, n_node) ] \ eye(n_node)
    else
        error("Element order not implemented.")
    end # end if
                
    
    return RefEl_Pk{n_order}(info_str, n_order, node, n_node, coeff)
    
end # end function


# -------------------------------------------------------------------------------------
# --------------------   Functions for RefEl{n_order}   --------------------
# -------------------------------------------------------------------------------------


# ---------------------------------------------
function shapeFun(ref_el :: RefEl_Pk{1}, p :: Array{Float64})

    size(p,1)!=2 ? error("Point array must be of size=2-by-n.") :

    value  = ref_el.coeff * [p ; ones(1, size(p,2))]
            
    return value
end


function shapeFun_grad(ref_el :: RefEl_Pk{1}, p :: Array{Float64})
    
    size(p,1)!=2 ? error("Point array must be of size=2-by-n.") :

    value = reshape([ref_el.coeff[i,1:2] for j=1:size(p,2) for i=1:ref_el.n_node], ref_el.n_node, size(p,2))
    
    return value
end
# ---------------------------------------------


# ---------------------------------------------
function shapeFun(ref_el :: RefEl_Pk{2}, p :: Array{Float64})
    
    size(p,1)!=2 ? error("Point array must be of size=2-by-n.") :

    value  = ref_el.coeff * [p.^2 ;  p[1,:]'.*p[1,:]' ; p ; ones(1,size(p,2))]
    
    return value
end


function shapeFun_grad(ref_el :: RefEl_Pk{2}, p :: Array{Float64})

    size(p,1)!=2 ? error("Point array must be of size=2-by-n.") :

    dx = ref_el.coeff * [2*p[1,:] ; zeros(1,size(p,2)) ; p[2,:]' ; ones(1,size(p,2)) ; zeros(1,size(p,2)) ; zeros(1,size(p,2))]
    dy = ref_el.coeff * [zeros(1,size(p,2)) ; 2*p[2,:] ; p[1,:]' ; zeros(1,size(p,2)) ; ones(1,size(p,2)) ; zeros(1,size(p,2))]
           
    value = reshape([[dx[i] ; dy[i]] for i in eachindex(dx)], ref_el.n_node, size(p,2))

    return value
end
# ---------------------------------------------


# ---------------------------------------------
function shapeFun(ref_el :: RefEl_Pk{3}, p :: Array{Float64})
    
    size(p,1)!=2 ? error("Point array must be of size=2-by-n.") :

    value  = ref_el.coeff * [ p.^3 ; p[1,:]'.^2.*p[2,:]' ; p[1,:]'.*p[2,:]'.^2 ; 
                                p[1,:].^2 ; p[2,:].^2 ; p[1,:]'.*p[2,:]' ; 
                                p ; 
                                ones(1, size(p,2)) ]
    
    return value
end


function shapeFun_grad(ref_el :: RefEl_Pk{3}, p :: Array{Float64})

    size(p,1)!=2 ? error("Point array must be of size=2-by-n.") :
    
    dx = ref_el.coeff * [ 3*p[1,:].^2 ; zeros(1, size(p,2)) ; 2*p[1,:]'.*p[2,:]' ; p[2,:]'.^2 ;
                            2*p[1,:] ; zeros(1,size(p,2)) ; p[2,:]' ; 
                            ones(1,size(p,2)) ; zeros(1,size(p,2)) ; 
                            zeros(1,size(p,2)) ]
    dy = ref_el.coeff * [ zeros(1, size(p,2)) ; 3*p[2,:].^2 ; p[1,:]'.^2 ; 2*p[1,:]'.*p[2,:]'
                            zeros(1,size(p,2)) ; 2*p[2,:] ; p[1,:]' ; 
                            zeros(1,size(p,2)) ; ones(1,size(p,2)) ; 
                            zeros(1,size(p,2)) ]
           
    value = reshape([[dx[i] ; dy[i]] for i in eachindex(dx)], ref_el.n_node, size(p,2))

    return value
end
# ---------------------------------------------

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
