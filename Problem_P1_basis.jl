type BasisFun <: AbstractBasisProblem
    
    info :: String
    type_info :: String
    
    T :: Float64

    coeff :: Array{Float64,2}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool

    problem_parent :: AbstractProblem
    
    function BasisFun(problem :: AbstractProblem, tri :: Geometry.Triangle)
        
        info = "Evolution of basis."
        type_info = problem.type_info

        T = problem.T

        coeff = Geometry.compute_P1_basis_coeff(tri)

        is_transient_diffusion = problem.is_transient_diffusion
        is_transient_velocity = problem.is_transient_velocity

        problem_parent = problem
        
        return new(info, type_info, T, coeff, is_transient_diffusion, is_transient_velocity, problem_parent)
    end # end constructor
end # end type



# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------


function diffusion(problem :: BasisFun, t :: Float64, x :: Array{Float64})
    """

    Delegate function calls to the original abstract problem.

    """
    
    return diffusion(problem.problem_parent, t, x)
end


# --------------------------------------------------------------------


function velocity(problem :: BasisFun,  t :: Float64, x :: Array{Float64})
    """

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

    """

    return velocity(problem.problem_parent, t, x)
end


# --------------------------------------------------------------------


function u_init(problem :: BasisFun, x :: Array{Float64,1})
                
    length(x)!=2 ? error(" Vector x must length=2.") :
   
    out = [x ; 1]'  * problem.coeff
                
    return out
end


function u_init(problem :: BasisFun, x :: Array{Float64,2})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :
    
    out = [x ones(size(x,1))] * problem.coeff
    
    return out
end

function u_init(problem :: BasisFun, x :: Array{Float64,3})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :
        
    out = Array{Float64, 3}(size(x,1), 3, size(x,3))

    for i=1:size(x,3)
        out[:,:,i] = u_init(problem, x[:,:,i])
    end
    
    return out
end
