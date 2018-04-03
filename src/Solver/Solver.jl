module Solver

import FEM, Mesh, Parameter, Problem, Quad, TimeIntegrator, Geometry, FiniteDiff, Reconstruction

using ProgressMeter


# -------   Core routines for workflow   -------
include("Solver_FEM.jl")
include("Solver_MsFEM_std.jl")
include("Solver_MsFEM_reconstruction.jl")


end