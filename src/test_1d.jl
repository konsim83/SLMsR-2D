import Mesh, Problem, FEM, Quad, TimeIntegrator, Plot

import DifferentialEquations


function traceback(point :: Array{Float64,2},
					T :: Float64, 
					velocity :: Function)
	dt = 1/10
	n_steps_f = 5

	tspan = (-T, -T+dt)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob, dtmax=dt/n_steps_f)

	return sol.u
end


# ----------------------------------------------------------------
m_coarse = Mesh.mesh_unit_square(2)
m_simplex = Mesh.mesh_unit_simplex()

# If this is not the case then the conformal MsFEM is just the
# standard FEM
m_simplex = Mesh.refine_rg(m_simplex, 4)

# Subdividing edges messes up boundary markers. We need to correct
# that.
ind_point_boundary = sort(unique(m_simplex.edge[:,m_simplex.edge_marker.!=0]))
m_simplex.point_marker[:] = zeros(Int, size(m_simplex.point_marker))
m_simplex.point_marker[ind_point_boundary] = ones(Int, size(ind_point_boundary))

mesh_collection = Mesh.TriMesh_collection(m_coarse, m_simplex)
# ----------------------------------------------------------------


problem = Problem.Gaussian_R_4_conserv(1.0, 0.05 , 30)


# ----------------------------------------------------------------
velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

	return hcat(-Problem.velocity(problem, -t, x)...)
end
# ----------------------------------------------------------------

T = 1.0

point_all = hcat([mesh.point for mesh in mesh_collection.mesh_f]...)
point_orig_all = traceback(point_all, T, velocity)

U = Problem.u_init(problem, mesh_collection.mesh.point)
u_before = Problem.u_init(problem, point_orig_all[1])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	i_mesh = 1 # mesh index
	j = 2 # Time
	seg = 3 # segment

	mesh_local = mesh_collection.mesh_f[i_mesh]

	
	cell_2d = circshift(mesh_local.segment[:,mesh_local.segment_marker.==seg],1)
	ind_edge = sort(unique(cell_2d))
	n = length(ind_edge)


	tri = Geometry.Triangle(mesh_collection.mesh.point[:, mesh_collection.mesh.cell[:,i_mesh]])
    problem_f = Problem.BasisFun(problem, tri)
    
	ind = (1:mesh_local.n_point) + (i_mesh-1)*mesh_local.n_point
	
	u_before_local = u_before[ind_edge]
	u0_local = Problem.u_init(problem_f, point_orig_all[1][:,ind])

	if seg==1
		ind_con = [1 ; 2 ; 1+n ; 2+n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[mesh_collection.mesh.cell[:,i_mesh]][1]
	    a_2 = U[mesh_collection.mesh.cell[:,i_mesh]][2]

		basis_left = u0_local[ind_edge,1]
		basis_right = u0_local[ind_edge,2]
	elseif seg==2
		ind_con = [1 ; 2 ; 1+n ; 2+n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[mesh_collection.mesh.cell[:,i_mesh]][2]
	    a_2 = U[mesh_collection.mesh.cell[:,i_mesh]][3]

		basis_left = u0_local[ind_edge,2]
		basis_right = u0_local[ind_edge,3]
	elseif seg==3
		ind_con = [1 ; 2 ; 1+n ; 2+n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[mesh_collection.mesh.cell[:,i_mesh]][1]
	    a_2 = U[mesh_collection.mesh.cell[:,i_mesh]][3]

		basis_left = u0_local[ind_edge,1]
		basis_right = u0_local[ind_edge,3]
	end


	k_edge = 0.01

	system_matrix = [ (a_1^2 + k_edge)*speye(n)   (a_1*a_2)*speye(n) ; 
                        (a_1*a_2)*speye(n)   (a_2^2 + k_edge)*speye(n) ]

    rhs = ( [k_edge*basis_left + a_1*u_before_local ;
            k_edge*basis_right + a_2*u_before_local]
            - system_matrix[:,ind_con]*val_con )


    basis_lr = zeros(2*n)
    basis_lr[ind_con] = val_con
    basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]








	
# 		# Connectivity
# 		cell_2d = circshift(mesh_local.segment[:,mesh_local.segment_marker.==seg],1)

# 		# Indices in mesh_local coordinates, not ordered
# 		ind_edge = sort(unique(cell_2d))
# 		n = length(ind_edge)

# 		# point list
# 		p_old = point_orig_all[end-(j-2)][:,ind][:,ind_edge]
# 		p_next = point_orig_all[end-(j-1)][:,ind][:,ind_edge]
# # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# for i=1:n
# 	cell_2d[cell_2d.==ind_edge[i]] = i
# end






