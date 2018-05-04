struct Gaussian_R_5_conserv <: AbstractPhysicalProblem
    
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

    psi :: Float64

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool

    conservative :: Bool

    k :: Int
    
end # end type


function Gaussian_R_5_conserv(T :: Float64, psi :: Float64, k :: Int)
        
    info_prob = "Evolution of symmetric Gaussian_R_5_conserv (classic velocity)."
    type_info = "ADE"
    file_name = "Gaussian_R_5_conserv"

    marker_dirichlet_edge = Array{Int}(0)
    marker_neumann_edge = Array{Int}(0)

    lambda_1 = 0.03
    lambda_2 = 0.03

    alpha = 0*pi/8
    rot = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]
    
    covariance_mat = rot * diagm([lambda_1 ; lambda_2]) * rot'
    covariance_mat_det = det(covariance_mat)
    covariance_mat_inv = covariance_mat \ eye(2)
    expectation = [1/2 ; 1/2]

    is_transient_diffusion = false
    is_transient_velocity = true

    conservative = true
    
    return Gaussian_R_5_conserv(info_prob, type_info, file_name,
                    T, 
                    marker_dirichlet_edge, marker_neumann_edge,
                    covariance_mat, covariance_mat_det, covariance_mat_inv, 
                    expectation,
                    psi, 
                    is_transient_diffusion, 
                    is_transient_velocity,
                    conservative,
                    k)
end # end constructor


# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    diffusion(problem :: Gaussian_R_5_conserv, t :: Float64, x :: Array{Float64,2})

    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian_R_5_conserv, t :: Float64, x :: Array{Float64,2})
    
    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    out = [[0.01*(1-0.9999*sin(2*pi*problem.k*x[1,i])) 0.0 ; 
            0.0 0.01*(1-0.9999*sin(2*pi*problem.k*x[2,i]))] for i=1:size(x,2)]
    
    return out
end


"""
    diffusion(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Array{Float64,2},1})
    
    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Array{Float64,2},1})
        
    out = [diffusion(problem, t, y) for y in x]
    
    return out
    
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------
function velocity_div(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    psi = problem.psi
    T = 1.

    V = (psi/T)*8*pi^2 * sin.(2*pi*(x[1,:]-t/T)) .* cos.(2*pi*(x[1,:]-t/T)) .* cos.(pi*(x[2,:]-1/2)) * ( 
            cos.(pi*(x[2,:]-1/2)) - sin.(pi*(x[2,:]-1/2))
          ) 
                

    out = [V[i,:] for i=1:size(x,2)]
    
    return out
end

function velocity_div(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    out = [velocity_div(problem, t, y) for y in x]

    return out
end





"""
    velocity(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Float64,2})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    psi = problem.psi
    T = 1.

    V = psi * hcat( (psi/T)*2*pi * (sin.(2*pi*(x[1,:]-t/T)).^2) .* cos.(pi*(x[2,:]-1/2)).^2, 
                    (psi/T)*4*pi * sin.(2*pi*(x[1,:]-t/T)) .* cos.(2*pi*(x[1,:]-t/T)) .* (cos.(pi*(x[2,:]-1/2)).^2)
                   )

    out = [V[i,:] for i=1:size(x,2)]
    
    return out
end


"""
    velocity(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: Gaussian_R_5_conserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    out = [velocity(problem, t, y) for y in x]

    return out
end



# --------------------------------------------------------------------
# --------------------------------------------------------------------

function u_init(problem :: Gaussian_R_5_conserv, x :: Array{Float64})
                
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :

    x = broadcast(+, -problem.expectation, x)
    x1 = broadcast(+, -[1/4 ; 0], x)
    x2 = broadcast(+, [1/4 ; 0], x)
    
    out  = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * (
                exp.( -1/2 * sum(x1.*(problem.covariance_mat_inv*x1),1) )
                + exp.( -1/2 * sum(x2.*(problem.covariance_mat_inv*x2),1) )
                ) / 2
    
    return vec(out)
end