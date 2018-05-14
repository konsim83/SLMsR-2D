module Problem


using Geometry


# -------   Type hierarchies   -------
abstract type AbstractProblem end
abstract type AbstractPhysicalProblem <: AbstractProblem end
abstract type AbstractBasisProblem <: AbstractProblem end

# ------------------------------------------------------------------------    
include("Problem_gaussian.jl")
include("Problem_gaussian_1.jl")
include("Problem_gaussian_2.jl")
include("Problem_gaussian_2a.jl")

include("Problem_gaussian_R_1.jl")

include("Problem_gaussian_R_2.jl")
include("Problem_gaussian_R_2_conserv.jl")

include("Problem_gaussian_R_3.jl")
include("Problem_gaussian_R_3_conserv.jl")

include("Problem_gaussian_R_4.jl")
include("Problem_gaussian_R_4_conserv.jl")

include("Problem_gaussian_R_5.jl")
include("Problem_gaussian_R_5_conserv.jl")

#include("Problem_wave.jl")
#include("Problem_rect_1.jl")
#include("Problem_hat_function_1.jl")

include("Problem_P1_basis.jl")

end # end module
