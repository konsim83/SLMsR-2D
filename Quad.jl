module Quad

# -------   Type hierarchies   -------
abstract type AbstractQuad end 

include("Quad_line.jl")
include("Quad_simplex.jl")

end # end module




