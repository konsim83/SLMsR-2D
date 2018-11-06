module Problem


using Geometry, Mesh, LinearAlgebra


# -------   Type hierarchies   -------
abstract type AbstractProblem end
abstract type AbstractPhysicalProblem <: AbstractProblem end
abstract type AbstractBasisProblem <: AbstractProblem end

# ------------------------------------------------------------------------    
include("Problem_GaussianSolenoidal.jl")

include("Problem_GaussianDivergent.jl")
include("Problem_GaussianDivergentConserv.jl")

include("Problem_GaussianRandomized.jl")

include("Problem_P1_basis.jl")

end # end module