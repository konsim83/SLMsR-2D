module PostProcess

import Mesh, FEM

using ProgressMeter

import TriangleMesh, NearestNeighbors


# --------------------------------------------
include("PostProcess_find_cell.jl")
# --------------------------------------------

# --------------------------------------------
include("PostProcess_evaluate_struct_FEM.jl")
include("PostProcess_evaluate_FEM.jl")

include("PostProcess_evaluate_struct_MsFEM.jl")
include("PostProcess_evaluate_MsFEM.jl")

include("PostProcess_evaluateGrad_FEM.jl")

include("PostProcess_map_solution.jl")
# --------------------------------------------

# --------------------------------------------
# Norms and graphs
include("PostProcess_error_L2.jl")
include("PostProcess_error_H1.jl")
# --------------------------------------------

end