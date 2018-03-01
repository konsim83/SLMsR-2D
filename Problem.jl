module Problem


using Geometry


# -------   Type hierarchies   -------
abstract type AbstractProblem end
abstract type AbstractPhysicalProblem <: AbstractProblem end
abstract type AbstractBasisProblem <: AbstractProblem end

# ------------------------------------------------------------------------    
include("Problem_gaussian.jl")
#include("Problem_gaussian_1.jl")
#include("Problem_gaussian_2.jl")
#include("Problem_gaussian_2a.jl")
#include("Problem_wave.jl")
#include("Problem_rect_1.jl")
#include("Problem_hat_function_1.jl")

include("Problem_P1_basis.jl")

end # end module
