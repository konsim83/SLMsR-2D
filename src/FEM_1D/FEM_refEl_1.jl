type RefEl_Lagrange_1
    info :: String
    order :: Int64

    node :: Array{Float64, 1}
    coeff :: Array{Float64, 2}

    eval :: Function
    eval_der :: Function

    function RefEl_Lagrange(order :: Int64)
        this = new()
        
        order != 1 ? error("Use order 1 only!") :

        this.info = string("Lagrange reference element of order ", string(order), " on intervall [0, 1].")
        this.order = order
        
        if order==0
            this.node = [0.5]
            this.coeff = [1]'
        else
            this.node = collect(linspace(0,1,order + 1))
            P = transpose(repmat(copy(this.node), 1, order + 1))
            for k = 1:(order+1)
                P[k,:] = P[k,:].^(order- (k-1))
            end
            this.coeff = eye(order + 1)/P
        end
        
        # ------------------------------------------------------------------    
        this.eval = function(x :: Array{Float64,1})
            
            if this.order==0
                val = ones(size(x))'
            else
                X = transpose(repmat(copy(x), 1, this.order + 1))
                for k = 1:(this.order+1)
                    X[k,:] = X[k,:].^(this.order - (k-1))
                end
                val = this.coeff*X
            end
            
            return val
        end # end function
        # ------------------------------------------------------------------        

        # ------------------------------------------------------------------    
        this.eval_der = function(x :: Array{Float64,1})

            X = transpose(repmat(copy(x), 1, this.order))
            for k = 1:this.order
                X[k,:] = (this.order - (k-1))*X[k,:].^(this.order -  k)
            end
            
            return this.coeff[:,1:end-1]*X
        end # end function
        # ------------------------------------------------------------------        

        return this
    end # end constructor
            
end # end type