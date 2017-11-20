type RefEl_Pk
    info :: String

    order :: Int64
    
    node :: Array{Float64, 2}
    n_node :: Int64

    coeff :: Array{Float64, 2}

    eval :: Function
    eval_grad :: Function
    
    function RefEl_Pk(order :: Int64)
        this = new()

        this.info = "Triangular Lagrange element of type P$order."

        if order==1
            # -------   P1   -------
            this.order = 1
            
            this.node = [0 0 ; 1 0 ; 0 1]
            this.n_node = 3

            # columns are coefficients of basis functions
            # phi(x) = ax + by + c
            this.coeff = [this.node ones(this.n_node)]\eye(3)
            
            this.eval = function(p :: Array{Float64, 2})
                value  = [p ones(size(p,1))] * this.coeff
                
                return value
            end
                
            this.eval_grad = function(p :: Array{Float64, 2})
                value  = cat( 3, repmat(this.coeff[1,:], 1, size(p,1))', repmat(this.coeff[2,:], 1, size(p,1))' )
                
                return value
            end
                
        elseif order==2
                # -------   P2   -------
                this.order = 2
                
                this.node = [0 0 ; 1 0 ; 0 1 ; 0.5 0 ; 0.5 0.5 ; 0 0.5]
                this.n_node = 6

                # columns are coefficients of basis functions
                # phi(x) = ax^2 + by^2 + cxy + dx + ey + f
                this.coeff = [this.node[:,1].^2   this.node[:,2].^2   this.node[:,1].*this.node[:,2]   this.node   ones(this.n_node)] \ eye(6)
            
                this.eval = function(p :: Array{Float64, 2})
                    value  = [p[:,1].^2   p[:,2].^2   p[:,1].*p[:,2]   p   ones(size(p,1))] * this.coeff
                
                    return value
                end
                
                this.eval_grad = function(p :: Array{Float64, 2})
                    dx = [2*p[:,1]   zeros(size(p,1))   p[:,2]   ones(size(p,1))   zeros(size(p,1))   zeros(size(p,1))] * this.coeff
                    dy = [zeros(size(p,1))   2*p[:,2]   p[:,1]   zeros(size(p,1))   ones(size(p,1))   zeros(size(p,1))] * this.coeff
                
                    return cat(3, dx, dy)
                end
                
        end # end if
                    
        
        return this
    end # end function    
end # end type
