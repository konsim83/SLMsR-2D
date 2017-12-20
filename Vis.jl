module Vis

using FEM, Mesh, WriteVTK

include("Vis_writeVTK_FEM.jl")

include("Vis_writeVTK_multiscale_solution_coarse.jl")
include("Vis_writeVTK_multiscale_basis.jl")

end # end module
