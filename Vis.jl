module Vis

using FEM, Mesh, WriteVTK, ProgressMeter

include("Vis_writeVTK_FEM.jl")

include("Vis_writeVTK_multiscale_solution_coarse.jl")
include("Vis_writeVTK_multiscale_solution_fine.jl")
include("Vis_writeVTK_multiscale_basis.jl")

end # end module
