module FEM

using Mesh, Parameter, Problem, Quad


abstract type AbstractRefEl end
abstract type AbstractDof end
abstract type AbstractSolution end


# -------   Reference elements   -------
include("FEM_RefEl.jl")


# -------   Dof types   -------
include("FEM_DoF_Pk.jl")
include("FEM_DoF_Pk_periodic.jl")


# -------   Solution types   -------
include("FEM_Solution.jl")


# -------   Assembly routines   -------
include("FEM_assemble.jl")

end # end module
