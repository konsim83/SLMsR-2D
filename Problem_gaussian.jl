type Gaussian <: Problem_top
    
    info :: String

    T :: Float64
    
    covariance_mat :: Array{Float64,2}
    covariance_mat_det :: Float64
    covariance_mat_sqrt_inv :: Array{Float64,2}
    expectation :: Array{Float64,1}
    
    diffusion:: Function
    velocity :: Function
    u_init :: Function

    function Gaussian(T :: Float64)

        this = new()
        
        this.info = "Evolution of symmetric Gaussian."

        this.T = T

        lambda_1 = 0.1
        lambda_2 = 0.1
        this.covariance_mat = diagm([lambda_1 ; lambda_2])
        this.covariance_mat_det = lambda_1 * lambda_2
        this.covariance_mat_sqrt_inv = sqrtm( this.covariance_mat \ eye(2))
        expectation = [1/2 ; 1/2]

        
        this.diffusion = function(t :: Float64, x :: Array{Float64,2})
            """
            Diffusion is represented by a positive 2-by-2 tensor.
            """

            size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
                
            a_11 = 0.1 * ones(1, 1, size(x,1))
            a_12 = zeros(1, 1, size(x,1))
            a_21 = zeros(1, 1, size(x,1))
            a_22 = 0.1 * ones(1, 1, size(x,1))
            
            out = [a_11 a_12 ; a_21 a_22]

            return out
        end

        this.velocity = function(t :: Float64, x :: Array{Float64,2})
            """
            Velocity is represented by a 2-vector. The solenoidal part
            can be represented by a stream function.
            """

            size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
            
            v_1 = zeros(size(x, 1))
            v_2 = zeros(size(x, 1))
            
            out = [v_1 v_2]
            
            return out
        end

        this.u_init = function(x :: Array{Float64,2})
                
            size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
                
            out = 1/sqrt(2*pi*this.covariance_mat_det) * exp( -1/2 * vec(sum(this.covariance_mat_sqrt_inv * x', 1)) )
            
            return out
        end

        return this
    end # end constructor
end # end type
