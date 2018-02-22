struct Evaluate_FEM

	n_time_point :: Int
	data :: Array{Array{Float64,2},1}
	local_function :: Array{Function,1}

end

# Outer constructor
function Evaluate_FEM(sol :: FEM.Solution_FEM, mesh :: Mesh.TriangleMesh.TriMesh)

	n_time_point = size(sol.u,2)

	data = Array{Array{Float64,2},1}(mesh.n_cell)
	local_function = Array{Function,1}(mesh.n_cell)

	for i in 1:mesh.n_cell
		data[i] = sol.u[mesh.cell[i,:],:]
		
		local_function[i] = function(bary :: Array{Float64,1})
			length(bary)!=3 ? error("Barycentirc coordinates vector must be of length 3.") :

			return data[i][1,:]*bary[1] + data[i][2,:]*bary[2] + data[i][3,:]*bary[3]
		end
	end

	return Evaluate_FEM(n_time_point, data, local_function)
end