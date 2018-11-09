function get_domain_points(x :: Array{Float64})
    
	# size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    # a = 0-1e-12
    # b = 1+1e-12

    # map to interval [-0, 1]
    x_unit = x .% 1
    
    x_unit[x_unit.<0] = 1 .+ x_unit[x_unit.<0]
    
    return x_unit
end


 
# Slow old verson, still used by FEM
"""
	find_cell(mesh :: Mesh.TriangleMesh.TriMesh, x :: Array{Float64,2}, warning = true :: Bool)

    For a given set of points 'x' find all the cells in 'mesh' that each point
    is in.

"""
function find_cell(mesh :: Mesh.TriangleMesh.TriMesh, 
					x :: Array{Float64}, 
					warning = false :: Bool)

	size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

	x = get_domain_points(x)

	if warning
		warn("Geometric predicates are not exact. The result could be an
		inacurate decision about a point being in a certain cell or not.")
	end

	P = [Mesh.get_point(mesh, mesh.cell) ; ones(1,3,mesh.n_cell)]
	X = [x ; ones(1,size(x,2))]

	# x_cell = zeros(Int, size(x,1))
	# x_bary_coord = zeros(size(x,1),3)
	x_cell = [Int[] for idx in 1:size(x,2)]
	x_bary_coord = [[Float64[]] for idx in 1:size(x,2)]
	for i in 1:size(x,2)
		pop!(x_bary_coord[i])
	end

	for i in 1:size(P,3)
		bary_coord = round.(inv(P[:,:,i]) * X,digits=9)

		in_triangle = (sum((bary_coord.>=0.0) .& (bary_coord.<=1.0),dims=1).==3)[:]

		# x_cell[in_triangle] = i
		# x_bary_coord[in_triangle,:] = bary_coord[in_triangle,:]
		map(xin->push!(xin, i), x_cell[in_triangle])

		ind = findall(in_triangle)
		for i in ind
			push!(x_bary_coord[i], bary_coord[:,i][:])
		end
	end

	return x_cell, x_bary_coord
end







function find_cell(mesh :: Mesh.TriangleMesh.TriMesh, 
					meshData :: Mesh.MeshData, 
					x :: Array{Float64}, 
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


"""
	find_cell(mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

    For a given set of points 'x' find all the cells in the mesh_collection,
    i.e., first the coarse cell and then for each coarse cell all fine cells.

"""
function find_cell(mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64})

	# size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

	# Find coarse cell
	# println("Searching points in coarse mesh...")
	x_cell, x_bary = PostProcess.find_cell(mesh_collection.mesh, mesh_collection.meshData, x)
	# println("...done.")

	# Find fine cell
	x_cell_f = [Array{Array{Int,1},1}(undef, 0) for i in 1:size(x,2)]
	x_bary_f = [Array{Array{Array{Float64,1},1},1}(undef, 0) for i in 1:size(x,2)]
	
	for i in 1:size(x,2)
		x_cell_loc, x_bary_loc = find_cell(mesh_collection.mesh_f[x_cell[i][1]], 
											mesh_collection.meshData_f[x_cell[i][1]],
											x[:,i])
		push!(x_cell_f[i], x_cell_loc[1])
		push!(x_bary_f[i], x_bary_loc[1])
	end

	return x_cell, x_bary, x_cell_f, x_bary_f
end


"""
	find_cell(mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

    For a given set of points 'x' find all the cells in the mesh_collection,
    i.e., first the coarse cell and then for each coarse cell all fine cells.

"""
function find_cell_old(mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64})

	# size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

	# Find coarse cell
	# println("Searching points in coarse mesh...")
	x_cell, x_bary = PostProcess.find_cell(mesh_collection.mesh, x)
	# println("...done.")

	# Find fine cell
	x_cell_f = [[Int[] for ixd_inner in 1:length(x_cell[idx])] for idx in 1:size(x,2)]
	x_bary_f = [[[Float64[]] for ixd_inner in 1:length(x_cell[idx])] for idx in 1:size(x,2)]

	# N = length(x_cell)
	# p = Progress(N, 0.01, "Searching points in coarse mesh...", 10)
	for idx in 1:size(x,1)
		for ixd_inner in 1:length(x_cell[idx])
			pop!(x_bary_f[idx][ixd_inner])
		end
		
		# next!(p)
	end

	# N = length(x_cell)
 #    p = Progress(N, 0.01, "Searching points in fine meshes...", 10)
	for idx in 1:length(x_cell)
		y = x[:,idx]
		for ixd_inner in 1:length(x_cell[idx])
			m = mesh_collection.mesh_f[x_cell[idx][ixd_inner]]

			x_cell_loc, x_bary_loc = find_cell(m, y, false)
			
			x_cell_f[idx][ixd_inner] = x_cell_loc[1]
			x_bary_f[idx][ixd_inner] = x_bary_loc[1]
		end

		# next!(p)
	end

	return x_cell, x_bary, x_cell_f, x_bary_f
end