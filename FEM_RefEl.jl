struct RefEl_Pk{n_order} <: AbstractRefEl
    info_str :: String

    n_order :: Int64
    
    node :: Array{Float64, 2}
    n_node :: Int64

    coeff :: Array{Float64, 2}
    
    function RefEl_Pk(::Val{n_order}) where {n_order}

        info_str = "Triangular Lagrange element of type P$n_order."

        if n_order==1
            # -------   P1   -------
            n_order = 1
            
            node = [0 0 ; 1 0 ; 0 1]
            n_node = 3

            # columns are coefficients of basis functions
            # phi(x) = ax + by + c
            coeff = [node ones(n_node)]\eye(3)

        elseif n_order==2
                # -------   P2   -------
                n_order = 2
                
                node = [0 0 ; 1 0 ; 0 1 ; 0.5 0 ; 0.5 0.5 ; 0 0.5]
                n_node = 6

                # columns are coefficients of basis functions
                # phi(x) = ax^2 + by^2 + cxy + dx + ey + f
                coeff = [node[:,1].^2   node[:,2].^2   node[:,1].*node[:,2]   node   ones(n_node)] \ eye(6)        

        end # end if
                    
        
        return new{Val{N}}(info_str, n_order, node, n_node, coeff)
        
    end # end function    
end # end type



# -------------------------------------------------------------------------------------
# --------------------   Functions for RefEl{n_order}   --------------------
# -------------------------------------------------------------------------------------


# ---------------------------------------------
function eval(ref_el :: RefEl_Pk{1}, p :: Array{Float64, 2})
    value  = [p ones(size(p,1))] * ref_el.coeff
            
    return value'
end

function eval(ref_el :: RefEl_Pk{1}, p :: Array{Float64, 1})
    
            
    return eval(ref_el, p')
end

function eval_grad(ref_el :: RefEl_Pk{1}, p :: Array{Float64, 1})
    
    
    return eval_grad(ref_el, p')
end

function eval_grad(ref_el :: RefEl_Pk{1}, p :: Array{Float64, 2})
    value  = cat( 3, repmat(ref_el.coeff[1,:], 1, size(p,1)), repmat(ref_el.coeff[2,:], 1, size(p,1)) )
    
    return value
end
# ---------------------------------------------


# ---------------------------------------------
function eval(ref_el :: RefEl_Pk{2}, p :: Array{Float64, 1})
    
    return eval(ref_el, p')
end

function eval(ref_el :: RefEl_Pk{2}, p :: Array{Float64, 2})
    value  = [p[:,1].^2   p[:,2].^2   p[:,1].*p[:,2]   p   ones(size(p,1))] * ref_el.coeff
    
    return value'
end

function eval_grad(ref_el :: RefEl_Pk{2}, p :: Array{Float64, 1})
                
    return eval_grad(ref_el, p')
end

function eval_grad(ref_el :: RefEl_Pk{2}, p :: Array{Float64, 2})
    dx = 2*p[:,1] 

    dx = [2*p[:,1]   zeros(size(p,1))   p[:,2]   ones(size(p,1))   zeros(size(p,1))   zeros(size(p,1))] * ref_el.coeff
    dy = [zeros(size(p,1))   2*p[:,2]   p[:,1]   zeros(size(p,1))   ones(size(p,1))   zeros(size(p,1))] * ref_el.coeff
                
    return cat(3, dx', dy')
end
# ---------------------------------------------

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
