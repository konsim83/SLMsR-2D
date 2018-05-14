module Mesh

using WriteVTK

import TriangleMesh, NearestNeighbors

# ---------------------------------------------------------------------------------------------

include("Mesh_TriMesh.jl")
include("Mesh_MeshData.jl")

include("Mesh_TriMesh_collection.jl")

# ---------------------------------------------------------------------------------------------

include("Mesh_std_domains.jl")
include("Mesh_refine.jl")

# ---------------------------------------------------------------------------------------------

include("Mesh_PeriodicityInfo.jl")

# ---------------------------------------------------------------------------------------------

include("Mesh_write.jl")

# ---------------------------------------------------------------------------------------------
end # end module Mesh
