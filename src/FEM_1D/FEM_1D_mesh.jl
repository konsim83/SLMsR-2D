struct Mesh_1D

	point :: Array{Float64,2}
	cell :: Array{Int,2}

	n_point :: Int
	n_cell :: Int

	function Mesh_1D(cell_edge :: Array{Int,2},
						point :: Array{Float64,2})

		ind_edge = unique(cell_edge)

		# Renumber 
		cell_renumbered = similar(cell_edge)
		for i=1:length(ind_edge)
		   cell_renumbered[cell_edge.==ind_edge[i]] = i
		end

		n_point = size(point,2)
		n_cell = size(cell_renumbered,2)

		return new(point, cell_renumbered, n_point, n_cell)
	end
end

# -------------------------------------------------------------------------------------
function map_ref_point(mesh :: Mesh_1D, x :: Array{Float64,1}, ind_c :: Array{Int64,1})

	P = [repmat(mesh.point[:,mesh.cell[1,ind]],1,length(x)) + [(mesh.point[1,mesh.cell[2,ind]]-mesh.point[1,mesh.cell[1,ind]])*transpose(x) ; 
											  (mesh.point[2,mesh.cell[2,ind]]-mesh.point[2,mesh.cell[1,ind]])*transpose(x)  ]
			for ind in ind_c ]

	return P
end
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
function map_ref_point_der(mesh :: Mesh_1D, x :: Array{Float64,1}, ind_c :: Array{Int64,1})

	P_der = sqrt.(sum((mesh.point[:,mesh.cell[2,ind_c]]-mesh.point[:,mesh.cell[1,ind_c]]).^2,1))

	return vec(P_der)
end
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
function get_cell_normal(mesh :: Mesh_1D, ind_c :: Array{Int64,1})

	rot = [0. 1.;-1. 0.]
	normal = [rot*(mesh.point[:,mesh.cell[2,i]]-mesh.point[:,mesh.cell[1,i]]) /
				sqrt(sum((mesh.point[:,mesh.cell[2,i]]-mesh.point[:,mesh.cell[1,i]]).^2)) for i in 1:mesh.n_cell]

	return normal
end
# -------------------------------------------------------------------------------------
