struct MeshData

	# Store the tree for the neighbor search
	tree :: NearestNeighbors.KDTree

	# For each point store the one ring as cell indices
	oneRingCells :: Array{Array{Int,1},1}

	# For each point store the inverse of all matrices (each cell of the one
	# ring has one matrix) to facilitate finding barycentric coordinates
	oneRingPointInv :: Array{Array{Array{Float64,2}}}

end

function MeshData(mesh :: TriangleMesh.TriMesh)

	tree = NearestNeighbors.KDTree(mesh.point)

	oneRingCells = []
	oneRingPointInv = []

	for i in 1:mesh.n_point
		
		if true
			oneRing = findall(vec(sum(mesh.cell.==i,dims=1)).!=0)

			push!(oneRingCells, oneRing)

		else
			# println("\n--------------------------------------")
			# println("Using two-ring for search structure...")
			# println("--------------------------------------\n")
			oneRing0 = findall(vec(sum(mesh.cell.==i,dims=1)).!=0)
			oneRing0 = unique(get_cell(mesh, oneRing0))
		
			oneRing = Array{Int,1}(0)
			for j=1:length(oneRing0)
				oneRing = [oneRing  ; findall(vec(sum(mesh.cell.==oneRing0[j],dims=1)).!=0)]
			end

			push!(oneRingCells, unique(oneRing))
		end

		P = [Mesh.get_point(mesh, Mesh.get_cell(mesh, oneRing)) ; ones(1,3,length(oneRing))]
		PInv = [inv(P[:,:,i]) for i in 1:size(P,3)]
		push!(oneRingPointInv, PInv)		
	end

	return MeshData(tree, oneRingCells, oneRingPointInv)
end