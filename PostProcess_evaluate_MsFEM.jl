function evaluate(E :: Evaluate_MsFEM, mesh :: Mesh.TriangleMesh.TriMesh, x :: Array{Float64,2})

	x_cell, x_bary = PostProcess.find_cell(mesh, x);

	u = [[zeros(E.n_time_point) for idx_inner in 1:length(x_cell[idx])] for idx in 1:size(x_cell,1)]

	for idx in 1:length(x_cell)
		for idx_inner in 1:length(x_cell[idx])
			f = E.local_function[x_cell[idx][idx_inner]]
			u[idx][idx_inner] = f(x_bary[idx][idx_inner])
		end
	end

	return u
end



function evaluate(E :: Evaluate_MsFEM, mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

	x_cell, x_bary, x_cell_f, x_bary_f = PostProcess.find_cell(mesh_collection, x);

	# u = [[zeros(E.n_time_point) for idx_inner in 1:length(x_cell[idx])] for idx in 1:size(x_cell,1)]

	# for idx in 1:length(x_cell)
	# 	for idx_inner in 1:length(x_cell[idx])
	# 		f = E.local_function[x_cell[idx][idx_inner]]
	# 		u[idx][idx_inner] = f(x_bary[idx][idx_inner])
	# 	end
	# end

	return u
end