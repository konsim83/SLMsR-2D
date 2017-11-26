type Gaussian <: Problem_top
    
    info :: String

    T :: Float64
    
    covariance_mat :: Array{Float64,2}
    covariance_mat_det :: Float64
    covariance_mat_sqrt_inv :: Array{Float64,2}
    expectation :: Array{Float64,1}

    
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

        return this
    end # end constructor
end # end type



# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------


function diffusion(problem :: Gaussian, t :: Float64, x :: Array{Float64,2})
    """

    Diffusion is represented by a positive 2-by-2 tensor.

    """
    
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2)
    
    out[:,1,1] = 0.1
    out[:,2,1] = 0.0
    out[:,1,2] = 0.0
    out[:,2,2] = 0.2 
    
    return out
end


function diffusion(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})
    """

    Diffusion is represented by a positive 2-by-2 tensor.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2, size(x,3))

    for i=1:size(x,3)
        out[:,:,:,i] = diffusion(problem, t, x[:,:,i])
    end
    
    return out
    
end


# --------------------------------------------------------------------


function velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,2})
    """

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    v_1 = zeros(size(x, 1))
    v_2 = zeros(size(x, 1))
            
    out = [v_1 v_2]
    
    return out
end


function velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})
    """

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64, 3}(size(x,1), 2, size(x,3))

    for i=1:size(x,3)
        out[:,:,i] = velocity(problem, t, x[:,:,i])
    end
    
    return out
    
end


# --------------------------------------------------------------------


function u_init(problem :: Gaussian, x :: Array{Float64,2})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :
                
    out = 1/sqrt(2*pi*problem.covariance_mat_det) * exp( -1/2 * vec(sum(problem.covariance_mat_sqrt_inv * x', 1)) )
    
    return out
end

function u_init(problem :: Gaussian, x :: Array{Float64,3})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :

    out = Array{Float64, 3}(size(x,1), 1, size(x,3))

    for i=1:size(x,3)
        out[:,1,i] = u_init(problem, x[:,:,i])
    end
    
    return out
end
