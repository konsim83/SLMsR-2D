"""
	struct Evaluate_MsFEM

	Structure contains a data of the solution in an ordered fashion and a set
	of local functions used to evaluate the solution at points given in
	barycentric coordinates.

"""
struct Evaluate_MsFEM

	n_time_point :: Int

	data :: Array{Array{Float64,2},1}
	local_function :: Array{Function,1}

	data_f :: Array{Array{Array{Float64,2},1},1}
	local_function_f :: Array{Array{Function,1},1}
end


"""
	Evaluate_MsFEM(sol :: FEM.Solution_MsFEM, mesh_collection :: Mesh.TriMesh_collection)

	Outer constructor for the struct 'Evaluate_MsFEM'.

"""
function Evaluate_MsFEM(sol :: FEM.Solution_MsFEM, mesh_collection :: Mesh.TriMesh_collection)

	mesh =  mesh_collection.mesh
	n_time_point = size(sol.u,2)

	data = Array{Array{Float64,2},1}(undef, mesh.n_cell)
	local_function = Array{Function,1}(undef, mesh.n_cell)

	for i in 1:mesh.n_cell
		data[i] = sol.u[mesh.cell[:,i],:]
		
		local_function[i] = function(bary :: Array{Float64,1})
			length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

			return data[i][1,:]*bary[1] + data[i][2,:]*bary[2] + data[i][3,:]*bary[3]
		end
	end


	data_f = [Array{Array{Float64,2},1}(undef, mesh_collection.n_elem_f[i]) for i in 1:mesh.n_cell]
	local_function_f = [Array{Function,1}(undef, mesh_collection.n_elem_f[i]) for i in 1:mesh.n_cell]

	for i in 1:mesh.n_cell
		coarse_weight = sol.u[mesh.cell[:,i],:]
		mesh_f = mesh_collection.mesh_f[i]
		for j in 1:mesh_f.n_cell			
			data_f[i][j] = ( coarse_weight[[1;1;1],:] .* sol.phi_1[i][mesh_f.cell[:,j],:]
							+ coarse_weight[[2;2;2],:] .* sol.phi_2[i][mesh_f.cell[:,j],:]
							+ coarse_weight[[3;3;3],:] .* sol.phi_3[i][mesh_f.cell[:,j],:] )

			local_function_f[i][j] = function(bary :: Array{Float64,1})
				length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

				return data_f[i][j][1,:]*bary[1] + data_f[i][j][2,:]*bary[2] + data_f[i][j][3,:]*bary[3]
			end
		end
	end

	return Evaluate_MsFEM(n_time_point, data, local_function, data_f, local_function_f)
end


function local_function(sol :: FEM.Solution_MsFEM,
						mesh :: Mesh.TriangleMesh.TriMesh,
						ind_cell :: Int,
						bary :: Array{Float64,1},
						k_time :: Int)

	length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

	data = sol.u[mesh.cell[:,ind_cell],k_time]

	return data[1]*bary[1] + data[2]*bary[2] + data[3]*bary[3]
end


function local_function_f(sol :: FEM.Solution_MsFEM,
							mesh_collection :: Mesh.TriMesh_collection,
							ind_cell :: Int,
							ind_cell_f :: Int,
							bary :: Array{Float64,1},
							k_time :: Int)

	length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

	coarse_weight = sol.u[mesh_collection.mesh.cell[:,ind_cell],k_time]
	mesh_f = mesh_collection.mesh_f[ind_cell]

	data_f = ( coarse_weight[[1;1;1]] .* sol.phi_1[ind_cell][mesh_f.cell[:,ind_cell_f],k_time]
				+ coarse_weight[[2;2;2]] .* sol.phi_2[ind_cell][mesh_f.cell[:,ind_cell_f],k_time]
				+ coarse_weight[[3;3;3]] .* sol.phi_3[ind_cell][mesh_f.cell[:,ind_cell_f],k_time] )

	return data_f[1]*bary[1] + data_f[2]*bary[2] + data_f[3]*bary[3]
end




# """
# 	Evaluate_MsFEM(sol :: FEM.Solution_MsFEM, mesh_collection :: Mesh.TriMesh_collection, k_time :: Int)

# 	Outer constructor for the struct 'Evaluate_MsFEM'.

# """
# function Evaluate_MsFEM(sol :: FEM.Solution_MsFEM, mesh_collection :: Mesh.TriMesh_collection, k_time :: Int)

# 	mesh =  mesh_collection.mesh
# 	n_time_point = size(sol.u,2)

# 	data = Array{Array{Float64,2},1}(mesh.n_cell)
# 	local_function = Array{Function,1}(mesh.n_cell)

# 	for i in 1:mesh.n_cell
# 		data[i] = sol.u[mesh.cell[:,i],k_time:k_time]
		
# 		local_function[i] = function(bary :: Array{Float64,1})
# 			length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

# 			return data[i][1,:]*bary[1] + data[i][2,:]*bary[2] + data[i][3,:]*bary[3]
# 		end
# 	end


# 	data_f = [Array{Array{Float64,2},1}(mesh_collection.n_elem_f[i]) for i in 1:mesh.n_cell]
# 	local_function_f = [Array{Function,1}(mesh_collection.n_elem_f[i]) for i in 1:mesh.n_cell]

# 	for i in 1:mesh.n_cell
# 		coarse_weight = sol.u[mesh.cell[:,i],k_time:k_time]
# 		mesh_f = mesh_collection.mesh_f[i]
# 		for j in 1:mesh_f.n_cell			
# 			data_f[i][j] = ( coarse_weight[[1;1;1],:] .* sol.phi_1[i][mesh_f.cell[:,j],k_time:k_time]
# 							+ coarse_weight[[2;2;2],:] .* sol.phi_2[i][mesh_f.cell[:,j],k_time:k_time]
# 							+ coarse_weight[[3;3;3],:] .* sol.phi_3[i][mesh_f.cell[:,j],k_time:k_time] )

# 			local_function_f[i][j] = function(bary :: Array{Float64,1})
# 				length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

# 				return data_f[i][j][1,:]*bary[1] + data_f[i][j][2,:]*bary[2] + data_f[i][j][3,:]*bary[3]
# 			end
# 		end
# 	end

# 	return Evaluate_MsFEM(n_time_point, data, local_function, data_f, local_function_f)
# end