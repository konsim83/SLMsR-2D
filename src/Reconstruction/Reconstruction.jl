module Reconstruction

import Mesh, Parameter, Problem, FEM, FEM_1D, Quad, TimeIntegrator, FiniteDiff, PostProcess

import DifferentialEquations

using ProgressMeter, SparseArrays, LinearAlgebra


function speye(n :: Int; typename=Float64 :: Type)

	return sparse(diagm(0=>ones(typename, n)))
end

# ---------------------
# Data Types

# ---------------------

include("Reconstruction_SemiLagrange.jl")

include("Reconstruction_SemiLagrange_evolve_edge.jl")

include("Reconstruction_SemiLagrange_reconstruct_edge_L2.jl")
include("Reconstruction_SemiLagrange_reconstruct_edge_L2_sumEqualOne.jl")

include("Reconstruction_SemiLagrange_reconstruct_L2_conformal.jl")
include("Reconstruction_SemiLagrange_reconstruct_L2_nonconformal.jl")
include("Reconstruction_SemiLagrange_reconstruct_H1_conformal.jl")
include("Reconstruction_SemiLagrange_reconstruct_H1_conformal_sumEqualOne.jl")

end