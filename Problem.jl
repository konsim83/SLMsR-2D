module Problem

# -------   Type hierarchies   -------
abstract type AbstractProblem end

# ------------------------------------------------------------------------    
include("Problem_gaussian.jl")
#include("Problem_wave.jl")
#include("Problem_rect_1.jl")
#include("Problem_hat_function_1.jl")

include("Problem_P1_basis.jl")

end # end module
