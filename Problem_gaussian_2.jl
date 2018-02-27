struct Gaussian_2 <: AbstractPhysicalProblem
    
    info_prob :: String
    type_info :: String
    
    T :: Float64

    k :: Int
    
    covariance_mat :: Array{Float64,2}
    covariance_mat_det :: Float64
    covariance_mat_inv :: Array{Float64,2}
    expectation :: Array{Float64,1}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool
    
    function Gaussian_2(T :: Float64, k :: Int)
        
        info_prob = "Evolution of symmetric Gaussian."
        type_info = "ADE"

        # --------------------------
        lambda_1 = 0.03
        lambda_2 = 0.03

        alpha = 0*pi/8
        rot = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]
        
        covariance_mat = rot * diagm([lambda_1 ; lambda_2]) * rot'
        covariance_mat_det = det(covariance_mat)
        covariance_mat_inv = covariance_mat \ eye(2)
        expectation = [1/2 ; 1/2]
        # --------------------------

        is_transient_diffusion = false
        is_transient_velocity = false
        
        return new(info_prob, type_info, T, k,
                    covariance_mat, covariance_mat_det, covariance_mat_inv, 
                    expectation, 
                    is_transient_diffusion, 
                    is_transient_velocity)
    end # end constructor
end # end type



# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    streamFun(problem :: Gaussian_2, t :: Float64, x :: Array{Float64,2})

    Stream function for velocity.
"""
function streamFun(problem :: Gaussian_2, t :: Float64, x :: Array{Float64,2})

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :

    # out = sin.(2*pi*problem.k*(x[:,1]-x[:,2]))
    out = sqrt.(sum(x.^2,2))

end

"""
    streamFunDer(problem :: Gaussian_2, t :: Float64, x :: Array{Float64,2})

    Stream function skew-derivative.
"""
function streamFunDer(problem :: Gaussian_2, t :: Float64, x :: Array{Float64,2})

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :

    # out = -problem.k * cat(2,cos.(2*pi*problem.k*(x[:,1]-x[:,2])), cos.(2*pi*problem.k*(x[:,1]-x[:,2])))
    out = -cat(2, -x[:,2]./sqrt.(sum(x.^2,2)), x[:,1]./sqrt.(sum(x.^2,2)))

    return out
end

# --------------------------------------------------------------------
# --------------------------------------------------------------------

"""

    diffusion(problem :: Gaussian_2, t :: Float64, x :: Array{Float64,2})
    
    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian_2, t :: Float64, x :: Array{Float64,2})
    
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2)

    #strFun = streamFun(problem, t, problem.k*x)

    out[:,1,1] = 0.01
    out[:,2,1] = 0
    out[:,1,2] = 0
    out[:,2,2] = 0.01
    
    return out
end


"""
    diffusion(problem :: Gaussian_2,  t :: Float64, x :: Array{Float64,3})

    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian_2,  t :: Float64, x :: Array{Float64,3})
    

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2, size(x,3))

    for i=1:size(x,3)
        out[:,:,:,i] = diffusion(problem, t, x[:,:,i])
    end
    
    return out
    
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------


"""
    velocity(problem :: Gaussian_2,  t :: Float64, x :: Array{Float64,2})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian_2,  t :: Float64, x :: Array{Float64,2})

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
    
    return streamFunDer(problem, t, x)
end


"""
    velocity(problem :: Gaussian_2,  t :: Float64, x :: Array{Float64,3})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian_2,  t :: Float64, x :: Array{Float64,3})

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :

    vel = Array{Float64}(size(x,1), 2, size(x,3))

    for i=1:size(x,3)
        vel[:,:,i] = velocity(problem, t, x[:,:,i])
    end 

    return vel
    
end



# --------------------------------------------------------------------
# --------------------------------------------------------------------


function u_init(problem :: Gaussian_2, x :: Array{Float64,1})
                
    length(x)!=2 ? error(" Vector x must length=2.") :

    x -= problem.expectation
                
    out = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * exp.( -1/2 * x'*problem.covariance_mat_inv * x )
    
    return out
end


function u_init(problem :: Gaussian_2, x :: Array{Float64,2})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :

    x = broadcast(+, -problem.expectation', x)
    
    out  = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * exp.( -1/2 * sum((x*problem.covariance_mat_inv).*x,2) )
    
    return out
end

function u_init(problem :: Gaussian_2, x :: Array{Float64,3})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :

    out = Array{Float64, 3}(size(x,1), 1, size(x,3))

    for i=1:size(x,3)
        out[:,1,i] = u_init(problem, x[:,:,i])
    end
    
    return out
end
