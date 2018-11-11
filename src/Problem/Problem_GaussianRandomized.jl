struct GaussianRandomized <: AbstractPhysicalProblem
    
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

    conservative :: Bool

    mesh_diff :: Mesh.TriangleMesh.TriMesh
    meshData_diff :: Mesh.MeshData
    random_coeff_diff :: Array{Float64,1}

    scale :: Float64
    
end # end type


function GaussianRandomized(T :: Float64; n_seg_per_edge = 25 :: Int, scale = 0.2 :: Float64)
        
    info_prob = "Evolution of symmetric Gaussian (divergent velocity)."
    type_info = "ADE"
    file_name = "GaussianRandomized"

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

    conservative = false

    # ------------------------------------
    # Evaluation mesh for diffusion data
    mesh_diff = Mesh.mesh_unit_square(n_seg_per_edge)
    meshData_diff = Mesh.MeshData(mesh_diff)
    r = rand(mesh_diff.n_cell) # cell based random data
    random_coeff_diff = 0.1*((r.-minimum(r)) / maximum(r.-minimum(r)) .+ 0.0001 )
    # ------------------------------------

    # Evaluation mesh for velocity data
    # mesh_vel = Mesh.mesh_unit_square(40)
    
    return GaussianRandomized(info_prob, type_info, file_name,
                                T, 
                                marker_dirichlet_edge, marker_neumann_edge,
                                covariance_mat, covariance_mat_det, covariance_mat_inv, 
                                expectation,
                                is_transient_diffusion, 
                                is_transient_velocity,
                                conservative,
                                mesh_diff,
                                meshData_diff,
                                random_coeff_diff,
                                scale)
end # end constructor


# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    diffusion(problem :: GaussianRandomized, t :: Float64, x :: Array{Float64,2})

    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: GaussianRandomized, t :: Float64, x :: Array{Float64,2})
    
    # size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    x_cell, x_bary_coord = find_cell(problem.mesh_diff, 
                                        problem.meshData_diff,
                                        x)

    out = [problem.random_coeff_diff[x_cell[idx][1]]*[1.0 0.0 ; 0.0 1.0] for idx in 1:length(x_cell)]
    
    return out
end


"""
    diffusion(problem :: GaussianRandomized,  t :: Float64, x :: Array{Array{Float64,2},1})
    
    Diffusion is represented by a positive 2-by-2 tensor.

"""
function diffusion(problem :: GaussianRandomized,  t :: Float64, x :: Array{Array{Float64,2},1})
        
    out = [diffusion(problem, t, y) for y in x]
    
    return out
    
end


# --------------------------------------------------------------------
# --------------------------------------------------------------------
"""
    velocity(problem :: GaussianRandomized,  t :: Float64, x :: Array{Float64,2})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: GaussianRandomized,  t :: Float64, x :: Array{Float64,2})

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    # out = [zeros(2) for i=1:size(x,2)]
    k1 = 1
    k2 = 1
    V = hcat(sin.(2*pi*k1*(x[1,:].-t)) .* cos.(2*pi*k2*(x[2,:])) *2*pi*k2,
                -cos.(2*pi*k1*(x[1,:].-t)) .* sin.(2*pi*k2*(x[2,:])) *2*pi*k1
            ) * problem.scale

    out = [V[i,:] for i=1:size(x,2)]
    
    return out
end


"""
    velocity(problem :: GaussianRandomized,  t :: Float64, x :: Array{Array{Float64,2},1})

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

"""
function velocity(problem :: GaussianRandomized,  t :: Float64, x :: Array{Array{Float64,2},1})

    out = [velocity(problem, t, y) for y in x]

    return out
end



# --------------------------------------------------------------------
# --------------------------------------------------------------------

function u_init(problem :: GaussianRandomized, x :: Array{Float64})
                
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