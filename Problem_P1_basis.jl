struct BasisFun <: AbstractBasisProblem
    
    info_prob :: String
    type_info :: String
    
    T :: Float64

    index_dirichlet_edge :: Array{Int}
    index_neumann_edge :: Array{Int}

    coeff :: Array{Float64,2}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool

    problem_parent :: AbstractProblem
    
    function BasisFun(problem :: AbstractProblem, tri :: Geometry.Triangle)
        
        info_prob = "Evolution of basis on triangular cell."
        type_info = problem.type_info

        T = problem.T

        index_dirichlet_edge = [1 ; 2 ; 3]
        index_neumann_edge = Array{Int}(0)

        coeff = Geometry.compute_P1_basis_coeff(tri)

        is_transient_diffusion = problem.is_transient_diffusion
        is_transient_velocity = problem.is_transient_velocity

        problem_parent = problem
        
        return new(info_prob, type_info, 
                    T,
                    index_dirichlet_edge, index_neumann_edge, 
                    coeff, 
                    is_transient_diffusion, is_transient_velocity,
                    problem_parent)
    end # end constructor
end # end type



# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------

"""
    diffusion(problem :: BasisFun, t :: Float64, x :: Array{Float64})

    Delegate function calls to the original abstract problem.

"""
function diffusion(problem :: BasisFun, t :: Float64, x :: Array{Float64})
    
    return diffusion(problem.problem_parent, t, x)
end


# --------------------------------------------------------------------

"""
    velocity(problem :: BasisFun,  t :: Float64, x :: Array{Float64})

    Delegate function calls to the original abstract problem.

"""
function velocity(problem :: BasisFun,  t :: Float64, x :: Array{Float64})

    return velocity(problem.problem_parent, t, x)
end


# --------------------------------------------------------------------


function u_init(problem :: BasisFun, x :: Array{Float64})
                
    size(x,1)!=2 ? error("Input points must be of size 2-by-n.") :
    
    out = problem.coeff * [x ; ones(1,size(x,2))]
    
    return out
end

function u_init(problem :: BasisFun, x :: Array{Array{Float64},1})
        
    out = [Array{Float64}(size(y,2)) for y in x]

    for y in x
        out[i] = u_init(problem, y)
    end
    
    return out
end


# --------------------------------------------------------------------


function u_dirichlet(problem :: BasisFun, x :: Array{Float64})

    size(x,1)!=2 ? error("Input points must be of size 2-by-n.") :
    
    out = problem.coeff * [x ; ones(1,size(x,2))]

end


# --------------------------------------------------------------------