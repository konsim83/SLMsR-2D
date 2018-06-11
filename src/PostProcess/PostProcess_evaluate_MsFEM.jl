"""
	evaluate(E :: Evaluate_MsFEM, mesh :: Mesh.TriangleMesh.TriMesh, x :: Array{Float64,2})

	Evaluate Multiscale solution at coarse level for a given set of input
	points.

"""
function evaluate(E :: Evaluate_MsFEM, mesh :: TriangleMesh.TriMesh, x :: Array{Float64,2})

	# Where are the new points found in the original mesh and what are their
	# barycentric coordinates?
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


"""
	evaluate(E :: Evaluate_MsFEM, mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

	Evaluate Multiscale solution at fine level for a given set of input
	points.

"""
function evaluate(E :: Evaluate_MsFEM, mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

	# Where are the new points found in the original mesh and what are their
	# barycentric coordinates?
	x_cell, x_bary, x_cell_f, x_bary_f = PostProcess.find_cell(mesh_collection, x);

	u_f = zeros(size(x,2), E.n_time_point)

	for idx in 1:length(x_cell)
		f = E.local_function_f[ x_cell[idx][1]][x_cell_f[idx][1][1] ]
		u_f[idx,:] = f( x_bary_f[idx][1][1] )
	end

	# u_f = [[[zeros(E.n_time_point) for idx_inner_2 in 1:length(x_cell_f[idx][idx_inner])] for idx_inner in 1:length(x_cell_f[idx])] for idx in 1:size(x_cell_f,1)]

	# for idx in 1:length(x_cell)
	# 	for idx_inner in 1:length(x_cell[idx])
	# 		for idx_inner_2 in 1:length(x_cell[idx][idx_inner])
	# 			f = E.local_function_f[ x_cell[idx][idx_inner]][x_cell_f[idx][idx_inner][idx_inner_2] ]
	# 			u_f[idx][idx_inner][idx_inner_2] = f( x_bary_f[idx][idx_inner][idx_inner_2] )
	# 		end
	# 	end
	# end

	return u_f
end


"""
	evaluate(sol :: FEM.Solution_MsFEM, 
					mesh_collection :: Mesh.TriMesh_collection, 
					x :: Array{Float64,2}, 
					k_time :: Int)

	Evaluate MsFEM solution at fine level directly at time index k_time.

"""
function evaluate(sol :: FEM.Solution_MsFEM, 
					mesh_collection :: Mesh.TriMesh_collection, 
					x :: Array{Float64,2}, 
					k_time :: Int)

	# Where are the new points found in the original mesh and what are their
	# barycentric coordinates?
	x_cell, x_bary, x_cell_f, x_bary_f = PostProcess.find_cell(mesh_collection, x)

	u_f = zeros(size(x,2))

	for idx in 1:size(x,2)
		# display(x_cell[idx])
		# display(x_cell_f[idx])
		u_f[idx] = local_function_f(sol, 
									mesh_collection, 
									x_cell[idx][1], 
									x_cell_f[idx][1][1], 
									x_bary_f[idx][1][1], 
									k_time)
	end

	return u_f
end


"""
	evaluate(sol :: FEM.Solution_MsFEM, mesh :: TriangleMesh.TriMesh, x :: Array{Float64,2}, k_time :: Int)

	Evaluate Multiscale solution at coarse level for a given set of input
	points.

"""
function evaluate(sol :: FEM.Solution_MsFEM, 
					mesh :: TriangleMesh.TriMesh, 
					x :: Array{Float64,2}, 
					k_time :: Int)

	# Where are the new points found in the original mesh and what are their
	# barycentric coordinates?
	x_cell, x_bary = PostProcess.find_cell(mesh, x);

	u = zeros(size(x,2))

	for idx in 1:size(x,2)
			u[idx] = local_function(sol, mesh, x_cell[idx][1], x_bary[idx][1], k_time)
	end

	return u
end