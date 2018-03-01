struct Gaussian <: AbstractPhysicalProblem
    
    info_prob :: String
    type_info :: String
    
    T :: Float64
    
    index_dirichlet_edge :: Array{Int}
    index_neumann_edge :: Array{Int}

    covariance_mat :: Array{Float64,2}
    covariance_mat_det :: Float64
    covariance_mat_inv :: Array{Float64,2}
    expectation :: Array{Float64,1}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool
    
end # end type


function Gaussian(T :: Float64)
        
    info_prob = "Evolution of symmetric Gaussian."
    type_info = "ADE"

    T = T

    index_dirichlet_edge = Array{Int}(0)
    index_neumann_edge = Array{Int}(0)

    lambda_1 = 0.05
    lambda_2 = 0.05

    alpha = 0*pi/8
    rot = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]
    
    covariance_mat = rot * diagm([lambda_1 ; lambda_2]) * rot'
    covariance_mat_det = det(covariance_mat)
    covariance_mat_inv = covariance_mat \ eye(2)
    expectation = [1/2 ; 1/2]

    is_transient_diffusion = false
    is_transient_velocity = false
    
    return new(info_prob, type_info, 
                T, 
                index_dirichlet_edge, index_neumann_edge,
                covariance_mat, covariance_mat_det, covariance_mat_inv, 
                expectation, 
                is_transient_diffusion, 
                is_transient_velocity)
end # end constructor


# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    diffusion(problem :: Gaussian, t :: Float64, x :: Array{Float64,2})

    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian, t :: Float64, x :: Array{Float64,2})
    
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2)

    
    out[:,1,1] = 0.01
    out[:,2,1] = 0.0
    out[:,1,2] = 0.0
    out[:,2,2] = 0.01
    
    return out
end


"""
    diffusion(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})
    
    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})
    

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
    velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,2})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,2})

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
    
    return zeros(size(x))
end


"""
    velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian,  t :: Float64, x :: Array{Float64,3})

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :

    return zeros(size(x))
    
end



# --------------------------------------------------------------------
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
