module Mesh

using PyPlot, WriteVTK

import TriangleMesh

# ---------------------------------------------------------------------------------------------

include("Mesh_TriMesh.jl")
include("Mesh_TriMesh_collection.jl")

# ---------------------------------------------------------------------------------------------

include("Mesh_std_domains.jl")

# ---------------------------------------------------------------------------------------------

include("Mesh_write.jl")
include("Mesh_plot.jl")
include("Mesh_maps.jl")

# ---------------------------------------------------------------------------------------------
end # end module Mesh
