module FEM

using SparseArrays, LinearAlgebra

import Mesh, Parameter, Problem, Quad


abstract type AbstractRefEl end

abstract type AbstractDof end
abstract type AbstractDofCollection end

abstract type AbstractSolution end


# -------   Reference elements   -------
include("FEM_RefEl.jl")


# -------   Dof types   -------
include("FEM_DoF_Pk.jl")
include("FEM_DoF_Pk_periodic.jl")
include("FEM_DoF_collection.jl")

include("FEM_DoF_maps.jl")


# -------   Solution types   -------
include("FEM_Solution.jl")


# -------   Assembly routines   -------
include("FEM_assemble.jl")

end # end module
