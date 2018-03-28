module Solver

import FEM, Mesh, Parameter, Problem, Quad, TimeIntegrator, Geometry, FiniteDiff, Reconstruction

using ProgressMeter

# -------   Core routines for workflow   -------
include("Solver_setup.jl")
include("Solver_pipeline.jl")

end
