type Gaussian <: AbstractPhysicalProblem
    
    info :: String
    type_info :: String
    
    T :: Float64
    
    covariance_mat :: Array{Float64,2}
    covariance_mat_det :: Float64
    covariance_mat_inv :: Array{Float64,2}
    expectation :: Array{Float64,1}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool
    
    function Gaussian(T :: Float64)

        this = new()
        
        this.info = "Evolution of symmetric Gaussian."
        this.type_info = "ADE"

        this.T = T

        lambda_1 = 0.05
        lambda_2 = 0.05

        alpha = 0*pi/8
        rot = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]
        
        this.covariance_mat = rot * diagm([lambda_1 ; lambda_2]) * rot'
        this.covariance_mat_det = det(this.covariance_mat)
        this.covariance_mat_inv = this.covariance_mat \ eye(2)
        this.expectation = [1/2 ; 1/2]

        is_transient_diffusion = false
        is_transient_velocity = false
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

    
    out[:,1,1] = 0.1001-0.1*sin.(2*pi*20*x[:,1])
    out[:,2,1] = 0.0
    out[:,1,2] = 0.0
    out[:,2,2] = 0.1001-0.1*cos.(2*pi*15*x[:,2])
    

    #=
    k = 25
    out[:,1,1] = 0.01
    out[:,2,1] = 1/k * sin.(k*(x[:,1]-x[:,2]))
    out[:,1,2] = -1/k * sin.(k*(x[:,1]-x[:,2]))
    out[:,2,2] = 0.01 
    =#
    
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
    
    return 0*ones(size(x)) #[cos.(k*(x[:,1]-x[:,2])) cos.(k*(x[:,1]-x[:,2]))]#0*ones(size(x))
end


function velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})
    """

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :

    #out = Array{Float64, 3}(size(x,1), 2, size(x,3))

    return 0*ones(size(x)) #[cos.(k*(x[:,1]-x[:,2])) cos.(k*(x[:,1]-x[:,2]))]
    
end


# --------------------------------------------------------------------

function u_init(problem :: Gaussian, x :: Array{Float64,1})
                
    length(x)!=2 ? error(" Vector x must length=2.") :

    x -= problem.expectation
                
    out = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * exp.( -1/2 * x'*problem.covariance_mat_inv * x )
    
    return out
end


function u_init(problem :: Gaussian, x :: Array{Float64,2})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :

    x = broadcast(+, -problem.expectation', x)
    
    out  = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * exp.( -1/2 * sum((x*problem.covariance_mat_inv).*x,2) )
    
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
