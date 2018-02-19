module PostProcess

using Mesh, FEM, ProgressMeter

export find_cell

function find_cell(mesh :: Mesh.TriangleMesh.TriMesh, x :: Array{Float64,2})

	P = [Mesh.get_point(mesh, mesh.cell) ones(3,1,mesh.n_cell)]
	X = [x ones(size(x,1))]

	x_cell = zeros(Int, size(x,1))
	x_bary_coord = zeros(size(x,1),3)

	for i in 1:size(P,3)
		bary_coord = X * inv(P[:,:,i])

		in_triangle = (sum((bary_coord.>=0.0) .& (bary_coord.<=1.0),2).==3)[:]

		x_cell[in_triangle] = i
		x_bary_coord[in_triangle,:] = bary_coord[in_triangle,:]
	end

	return x_cell, x_bary_coord
end


function find_cell(mesh_collection :: Mesh.TriMesh_collection, x :: Array{Float64,2})

	error("Not implemented yet.")

	P = [Mesh.get_point(mesh, mesh.cell) ones(3,1,mesh.n_cell)]
	X = [x ones(size(x,1))]

	x_cell = zeros(Int, size(x,1))
	x_bary_coord = zeros(size(x,1),3)

	for i in 1:size(P,3)
		bary_coord = X * inv(P[:,:,i])

		in_triangle = (sum((bary_coord.>=0.0) .& (bary_coord.<=1.0),2).==3)[:]

		x_cell[in_triangle] = i
		x_bary_coord[in_triangle,:] = bary_coord[in_triangle,:]
	end

	return x_cell, x_bary_coord
end


end