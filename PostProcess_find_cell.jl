function find_cell(mesh :: Mesh.TriangleMesh.TriMesh, x :: Array{Float64,2}, warning = true :: Bool)

	if warning
		warn("Geometric predicates are not exact. The result could be an
		inacurate decision about a point being in a certain cell or not.")
	end

	P = [Mesh.get_point(mesh, mesh.cell) ones(3,1,mesh.n_cell)]
	X = [x ones(size(x,1))]

	# x_cell = zeros(Int, size(x,1))
	# x_bary_coord = zeros(size(x,1),3)
	x_cell = [Int[] for idx in 1:size(x,1)]
	x_bary_coord = [[Float64[]] for idx in 1:size(x,1)]
	for i in 1:size(x,1)
		pop!(x_bary_coord[i])
	end

	for i in 1:size(P,3)
		bary_coord = round.(X * inv(P[:,:,i]),12)

		in_triangle = (sum((bary_coord.>=0.0) .& (bary_coord.<=1.0),2).==3)[:]

		# x_cell[in_triangle] = i
		# x_bary_coord[in_triangle,:] = bary_coord[in_triangle,:]
		map(xin->push!(xin, i), x_cell[in_triangle])
		ind = find(in_triangle)
		for i in ind
			push!(x_bary_coord[i], bary_coord[i,:][:])
		end
	end

	return x_cell, x_bary_coord
end


function find_cell(mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

	# Find coarse cell
	x_cell, x_bary = PostProcess.find_cell(mesh_collection.mesh, x)

	# Find fine cell
	x_cell_f = [[Int[] for ixd_inner in 1:length(x_cell[idx])] for idx in 1:size(x,1)]
	x_bary_f = [[[Float64[]] for ixd_inner in 1:length(x_cell[idx])] for idx in 1:size(x,1)]
	for idx in 1:size(x,1)
		for ixd_inner in 1:length(x_cell[idx])
			pop!(x_bary_f[idx][ixd_inner])
		end
	end

	for idx in 1:length(x_cell)
		y = convert(Array{Float64,2}, x[idx,:]')
		for ixd_inner in 1:length(x_cell[idx])
			m = mesh_collection.mesh_f[x_cell[idx][ixd_inner]]

			x_cell_loc, x_bary_loc = find_cell(m, y, false)
			
			x_cell_f[idx][ixd_inner] = x_cell_loc[1]
			x_bary_f[idx][ixd_inner] = x_bary_loc[1]
		end
	end

	return x_cell, x_bary, x_cell_f, x_bary_f
end