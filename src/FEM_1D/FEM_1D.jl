module FEM_1D

import Quad, Problem

using LinearAlgebra, SparseArrays

include("FEM_1D_mesh.jl")
include("FEM_1D_refEl.jl")
include("FEM_1D_dof.jl")


include("FEM_1D_assemble.jl")

end