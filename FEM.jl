module FEM

using Mesh, Parameter

# -------   Reference elements   -------
include("FEM_RefEl_abstract.jl")
include("FEM_RefEl.jl")

# -------   Dof types   -------
include("FEM_DoF_abstract.jl")
include("FEM_DoF_Pk.jl")
include("FEM_DoF_Pk_periodic.jl")

# -------   Solution types   -------
include("FEM_Solution_abstract.jl")
include("FEM_Solution.jl")

end # end module
