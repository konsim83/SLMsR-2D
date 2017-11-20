module FEM

using Mesh

# -------   Reference elements   -------
include("FEM_RefEl.jl")

# -------   DoF types   -------
include("FEM_DoF_Pk.jl")
include("FEM_DoF_Pk_periodic.jl")


end # end module
