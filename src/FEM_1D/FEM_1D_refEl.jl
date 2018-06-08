struct RefEl_Lagrange_1
    
    order :: Int64
    node :: Array{Float64, 1}
    coeff :: Array{Float64, 2}

    function RefEl_Lagrange_1()
        
        order = 1
        node = [0.0 ; 1.0]        
        
        coeff = [-1.0 1.0 ; 
                  1.0 0.0  ]
        
        return new(order, node, coeff)
    end # end constructor
            
end # end type


# ------------------------------------------------------------------
function shapeFun(r :: RefEl_Lagrange_1, x :: Array{Float64,1})
    
    return r.coeff * transpose([x ones(length(x))])
end # end function
# ------------------------------------------------------------------


# ------------------------------------------------------------------
function shapeFunDer(r :: RefEl_Lagrange_1, x :: Array{Float64,1})
    
    return transpose([ones(length(x)) ; -ones(length(x))])
end # end function
# ------------------------------------------------------------------