module FEM

using Mesh, Parameter, Problem, Quad, Time_integrator


abstract AbstractRefEl
abstract AbstractDof
abstract AbstractSolution


# -------   Reference elements   -------
include("FEM_RefEl.jl")


# -------   Dof types   -------
include("FEM_DoF_Pk.jl")
include("FEM_DoF_Pk_periodic.jl")


# -------   Solution types   -------
include("FEM_Solution.jl")


# -------   Assembly routines   -------
include("FEM_assemble.jl")


# -------   Core routines for workflow   -------
include("FEM_core.jl")
include("FEM_Solver_pipeline.jl")

end # end module
