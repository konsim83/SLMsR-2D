module Solver

using FEM, Mesh, Parameter, Problem, Quad, Time_integrator, ProgressMeter


# -------   Core routines for workflow   -------
include("Solver_setup.jl")
include("Solver_pipeline.jl")

end
