module Vis

using FEM, Mesh, WriteVTK, ProgressMeter, LinearAlgebra, DelimitedFiles

include("Vis_writeVTK_FEM.jl")

include("Vis_writeVTK_multiscale_solution_coarse.jl")
include("Vis_writeVTK_multiscale_solution_fine.jl")
include("Vis_writeVTK_multiscale_basis.jl")

include("Vis_plotError.jl")
include("Vis_writeError.jl")

end # end module
