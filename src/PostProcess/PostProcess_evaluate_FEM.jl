function evaluate(E :: Evaluate_FEM, mesh :: TriangleMesh.TriMesh, x :: Array{Float64,2})

	x_cell, x_bary = PostProcess.find_cell(mesh, x)

	u = [[zeros(E.n_time_point) for idx_inner in 1:length(x_cell[idx])] for idx in 1:size(x_cell,1)]

	for idx in 1:length(x_cell)
		for idx_inner in 1:length(x_cell[idx])
			f = E.local_function[x_cell[idx][idx_inner]]
			
			u[idx][idx_inner] = f(x_bary[idx][idx_inner])
		end
	end

	return u
end


function evaluate(sol :: FEM.Solution_FEM, mesh :: TriangleMesh.TriMesh, x :: Array{Float64,2}, k_time :: Int)

	x_cell, x_bary = PostProcess.find_cell(mesh, x)

	u = zeros(size(x_cell,1))

	for idx in 1:length(x_cell)			
		u[idx] = local_function(sol, mesh, x_cell[idx][1], x_bary[idx][1], k_time)
	end

	return u
end