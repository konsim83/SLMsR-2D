module Problem

using Geometry, LinearAlgebra

import Mesh, NearestNeighbors


# -------   Type hierarchies   -------
abstract type AbstractProblem end
abstract type AbstractPhysicalProblem <: AbstractProblem end
abstract type AbstractBasisProblem <: AbstractProblem end


function get_domain_points(x :: Array{Float64})
    
	# size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    # map to interval [-0, 1]
    x_unit = x .% 1
    
    x_unit[x_unit.<0] = 1 .+ x_unit[x_unit.<0]
    
    return x_unit
end


function find_cell(mesh :: Mesh.TriangleMesh.TriMesh, 
                    meshData :: Mesh.MeshData, 
                    x :: Array{Float64};
                    warning = false :: Bool)

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    x = get_domain_points(x)
    x = x[:,:] # make it and 2d array

    if warning
        warn("Geometric predicates are not exact. The result could be an
        inacurate decision about a point being in a certain cell or not.")
    end

    # Find indices of nearest neighbor in mesh
    n_knn = 15
    idx, dist = NearestNeighbors.knn(meshData.tree, x, n_knn, true)


    x_bary_coord = Array{Array{Array{Float64,1},1},1}(undef, 0)
    x_cell = Array{Array{Int,1},1}(undef, 0)
    for i in 1:size(x,2)
        bary_coord = []
        cell = []
        counter = 1
        while isempty(cell) && counter<=n_knn
            oneRing = meshData.oneRingCells[idx[i][counter]...]
            PInv = meshData.oneRingPointInv[idx[i][counter]...]
            for j in 1:length(oneRing)
                b_coord = round.(PInv[j] * [x[:,i];1],digits=9)

                condition = all(0.0.<=b_coord.<=1.0)
                if condition
                    push!(bary_coord, b_coord)
                    push!(cell, oneRing[j])
                end
            end
            counter += 1
        end
        
        push!(x_bary_coord, bary_coord)
        push!(x_cell, cell)
    end

    return x_cell, x_bary_coord
end

# ------------------------------------------------------------------------    

include("Problem_GaussianSolenoidal.jl")

include("Problem_GaussianDivergent.jl")
include("Problem_GaussianDivergentConservative.jl")

include("Problem_GaussianRandomized.jl")

include("Problem_P1_basis.jl")

end # end module