module RefEl

type RefEl_Pk
    info :: String

    node :: Array{Float64, 2}
    n_node :: Int64

    coeff :: Array{Float64, 2}

    eval_basis :: Function
    eval_basis_der :: Function
    
    function RefEl_Pk()
        this = new()

        this.info = "linear 3-node triangular Lagrange element"

        this.node = [0 0 ; 1 0 ; 0 1]
        this.n_node = 3

        # columns are coefficients of basis functions
        this.coeff = [this.node ones(this.n_node)]\eye(3)

        this.eval_basis = function(p :: Array{Float64, 2})
            value  = [p ones(size(p,1))] * this.coeff
            
            return value
        end

        this.eval_basis_der = function(p :: Array{Float64, 2})
            value  = cat( 3, repmat(this.coeff[1,:], 1, size(p,1))', repmat(this.coeff[2,:], 1, size(p,1))' )
            
            return value
        end
        
        return this
    end # end function    
end # end type

end # end module
