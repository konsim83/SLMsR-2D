struct Gaussian_2a <: AbstractPhysicalProblem
    
    info_prob :: String
    type_info :: String
    file_name :: String
    
    T :: Float64
    
    marker_dirichlet_edge :: Array{Int}
    marker_neumann_edge :: Array{Int}

    covariance_mat :: Array{Float64,2}
    covariance_mat_det :: Float64
    covariance_mat_inv :: Array{Float64,2}
    expectation :: Array{Float64,1}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool

    k :: Int
    
end # end type


function Gaussian_2a(T :: Float64, k :: Int)
        
    info_prob = "Evolution of symmetric Gaussian_2a (stream function)."
    type_info = "ADE"
    file_name = "Gaussian_2a"

    T = T

    marker_dirichlet_edge = Array{Int}(0)
    marker_neumann_edge = Array{Int}(0)

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
    

    return Gaussian_2a(info_prob, type_info, file_name,
                    T, 
                    marker_dirichlet_edge, marker_neumann_edge,
                    covariance_mat, covariance_mat_det, covariance_mat_inv, 
                    expectation, 
                    is_transient_diffusion, 
                    is_transient_velocity,
                    k)
end # end constructor


# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    streamFun(problem :: Gaussian_2a, t :: Float64, x :: Array{Float64,2})

    Stream function for velocity.
"""
function streamFun(problem :: Gaussian_2a, t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    out = sin.(2*pi*problem.k*(x[1,:]-x[2,:]))

    return out
end

"""
    streamFunDer(problem :: Gaussian_2a, t :: Float64, x :: Array{Float64,2})

    Stream function skew-derivative.
"""
function streamFunDer(problem :: Gaussian_2a, t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    V = problem.k*[2*pi*cos.(2*pi*problem.k*(x[1,:]-x[2,:])) 2*pi*cos.(2*pi*problem.k*(x[1,:]-x[2,:]))]
    out = [V[i,:] for i=1:size(x,2)]

    return out
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------


"""
    diffusion(problem :: Gaussian_2a, t :: Float64, x :: Array{Float64,2})

    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian_2a, t :: Float64, x :: Array{Float64,2})
    
    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    strFun = streamFun(problem, t, x)
    out = [[0.01 -strFun[i] ; strFun[i] 0.01] for i=1:size(x,2)]
    
    # out = [[0.01 0.0 ; 0.0 0.01] for i=1:size(x,2)]

    return out
end


"""
    diffusion(problem :: Gaussian_2a,  t :: Float64, x :: Array{Array{Float64,2},1})
    
    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian_2a,  t :: Float64, x :: Array{Array{Float64,2},1})
        
    out = [diffusion(problem, t, y) for y in x]
    
    return out
    
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------


"""
    velocity(problem :: Gaussian_2a,  t :: Float64, x :: Array{Float64,2})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian_2a,  t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    # out = streamFunDer(problem, t, x)

    out = [[0.0 ; 0.0] for i=1:size(x,2)]
    
    return out
end


"""
    velocity(problem :: Gaussian_2a,  t :: Float64, x :: Array{Array{Float64,2},1})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian_2a,  t :: Float64, x :: Array{Array{Float64,2},1})

    out = [velocity(problem, t, y) for y in x]

    return out
end



# --------------------------------------------------------------------
# --------------------------------------------------------------------

function u_init(problem :: Gaussian_2a, x :: Array{Float64})
                
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :

    x = broadcast(+, -problem.expectation, x)
    
    out  = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * exp.( -1/2 * sum(x.*(problem.covariance_mat_inv*x),1) )
    
    return vec(out)
end