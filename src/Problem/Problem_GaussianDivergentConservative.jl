struct GaussianDivergentConserv <: AbstractPhysicalProblem
    
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
    k1 :: Int
    k2 :: Int

end # end type


function GaussianDivergentConserv(T :: Float64;
                            k = 30 :: Int,
                            k1 = 1 :: Int,
                            k2 = 1 :: Int)
        
    info_prob = "Evolution of symmetric Gaussian (divergent velocity, conservative)."
    type_info = "ADE"
    file_name = "GaussianDivergentConserv"

    marker_dirichlet_edge = Array{Int}(undef, 0)
    marker_neumann_edge = Array{Int}(undef, 0)

    lambda_1 = 0.03
    lambda_2 = 0.03

    alpha = 0*pi/8
    rot = [cos(alpha) sin(alpha) ; -sin(alpha) cos(alpha)]
    
    covariance_mat = rot * diagm(0=>[lambda_1 ; lambda_2]) * rot'
    covariance_mat_det = det(covariance_mat)
    covariance_mat_inv = covariance_mat \ I
    expectation = [1/2 ; 1/2]

    is_transient_diffusion = false
    is_transient_velocity = true

    conservative = true
    
    return GaussianDivergentConserv(info_prob, type_info, file_name,
                    T, 
                    marker_dirichlet_edge, marker_neumann_edge,
                    covariance_mat, covariance_mat_det, covariance_mat_inv, 
                    expectation,
                    is_transient_diffusion, 
                    is_transient_velocity,
                    conservative,
                    k, k1, k2)
end # end constructor


# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    diffusion(problem :: GaussianDivergentConserv, t :: Float64, x :: Array{Float64,2})

    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: GaussianDivergentConserv, t :: Float64, x :: Array{Float64,2})
    
    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    out = [[0.01*(1-0.9999*sin(2*pi*problem.k*x[1,i])) 0.0 ; 
            0.0 0.01*(1-0.9999*sin(2*pi*problem.k*x[2,i]))] for i=1:size(x,2)]
    
    return out
end


"""
    diffusion(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Array{Float64,2},1})
    
    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Array{Float64,2},1})
        
    out = [diffusion(problem, t, y) for y in x]
    
    return out
    
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------


"""
    velocity(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Float64,2})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    k1 = 1
    k2 = 1
    V = hcat(-cos.(2*pi*k1*(x[1,:].-t)) .* sin.(2*pi*k2*(x[2,:])) *2*pi*k2,
                -cos.(2*pi*k1*(x[1,:].-t)) .* cos.(2*pi*k2*(x[2,:])) .* sin.(2*pi*2*x[1,:])
            )

    rotation = [cos(2*pi*t)   sin(2*pi*t) ; 
                -sin(2*pi*t)   cos(2*pi*t)]
    out = [rotation*V[i,:] for i=1:size(x,2)]
    
    return out
end


"""
    velocity(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    out = [velocity(problem, t, y) for y in x]

    return out
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------
"""
    reaction(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Float64,2})

    Divergence of velocity field.

"""
function reaction(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    psi = problem.psi

    out = psi*8*pi^2 * sin.(2*pi*(x[1,:].-t)) .* cos.(2*pi*(x[1,:].-t)) .* cos.(pi*(x[2,:].-1/2)) .* (cos.(pi*(x[2,:].-1/2)) - sin.(pi*(x[2,:].-1/2)))
    
    return out
end


"""
    reaction(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function reaction(problem :: GaussianDivergentConserv,  t :: Float64, x :: Array{Array{Float64,2},1})

    out = [reaction(problem, t, y) for y in x]

    return out
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------


function u_init(problem :: GaussianDivergentConserv, x :: Array{Float64})
                
    size(x,1)!=2 ? error(" List of vectors x must be of size 2-by-n.") :

    x = broadcast(+, -problem.expectation, x)
    x1 = broadcast(+, -[1/4 ; 0], x)
    x2 = broadcast(+, [1/4 ; 0], x)
    
    out  = 1/sqrt((2*pi)^2*problem.covariance_mat_det) * (
                exp.( -1/2 * sum(x1.*(problem.covariance_mat_inv*x1),dims=1) )
                + exp.( -1/2 * sum(x2.*(problem.covariance_mat_inv*x2),dims=1) )
                ) / 2
    
    return vec(out)
end