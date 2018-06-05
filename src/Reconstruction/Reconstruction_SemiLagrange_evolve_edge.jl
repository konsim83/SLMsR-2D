function evolve_edge!(mesh_local  :: Mesh.TriangleMesh.TriMesh,
						uBasisOpt :: Array{Float64,2},
						point_orig :: Array{Array{Float64,2},1},
						problem_local :: Problem.AbstractBasisProblem,
						dt_f :: Float64, n_steps_f :: Int,
 						seg :: Int)

	# Connectivity
	cell_2d = circshift(mesh_local.segment[:,mesh_local.segment_marker.==seg],1)

	# Indices in mesh_local coordinates, not ordered
	ind_edge = sort(unique(cell_2d))
	n = length(ind_edge)

	# point list
	p_old = point_orig[end-(j-2)][:,ind][:,ind_edge]
	p_next = point_orig[end-(j-1)][:,ind][:,ind_edge]

	error("Stop here...")
	u_basis_tmp[:,:,j] = predict_boundary_value()

end