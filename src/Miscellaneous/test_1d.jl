import Mesh, Problem, FEM, Quad, TimeIntegrator, Plot

import DifferentialEquations

using PyPlot

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

close()

# ----------------------------------------------------------------
m_coarse = Mesh.mesh_unit_square(2)
m_simplex = Mesh.mesh_unit_simplex(10, 0.01)

#-----------------------------------------------------------------
# If this is not the case then the conformal MsFEM is just the
# standard FEM
# m_simplex = Mesh.refine_rg(m_simplex, 4)

# Subdividing edges messes up boundary markers. We need to correct
# that.
# ind_point_boundary = sort(unique(m_simplex.edge[:,m_simplex.edge_marker.!=0]))
# m_simplex.point_marker[:] = zeros(Int, size(m_simplex.point_marker))
# m_simplex.point_marker[ind_point_boundary] = ones(Int, size(ind_point_boundary))
#-----------------------------------------------------------------

mesh_collection = Mesh.TriMesh_collection(m_coarse, m_simplex)
# ----------------------------------------------------------------


problem = Problem.Gaussian_R_4_conserv(1.0, 0.05 , 30)


# ----------------------------------------------------------------
velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

	return hcat(-Problem.velocity(problem, -t, x)...)
end
# ----------------------------------------------------------------

T = 1.0
i_mesh = i_mesh+1 # mesh index
j = 2 # Time
seg = 1 # segment

mesh = mesh_collection.mesh_f[i_mesh]

tri = Geometry.Triangle(mesh_collection.mesh.point[:, mesh_collection.mesh.cell[:,i_mesh]])
problem_f = Problem.BasisFun(problem, tri)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	point_orig = traceback(mesh_collection.mesh_f[i_mesh].point, T, velocity)
	pOrig = point_orig[1]
	u0 = Problem.u_init(problem_f, pOrig)

	U = Problem.u_init(problem, mesh_collection.mesh.point)
	
	uOrig = Problem.u_init(problem, pOrig)
	
	cell = circshift(mesh.segment[:,mesh.segment_marker.==seg],1)
	ind_edge = unique(cell)
	n = length(ind_edge)	

	pOrigEdge = pOrig[:,ind_edge]
	uOrigEdge = uOrig[ind_edge]
	
	

	if seg==1
		ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[mesh_collection.mesh.cell[:,i_mesh]][1]
	    a_2 = U[mesh_collection.mesh.cell[:,i_mesh]][2]

		basis_left = u0[ind_edge,1]
		basis_right = u0[ind_edge,2]
	elseif seg==2
		ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[mesh_collection.mesh.cell[:,i_mesh]][2]
	    a_2 = U[mesh_collection.mesh.cell[:,i_mesh]][3]

		basis_left = u0[ind_edge,2]
		basis_right = u0[ind_edge,3]
	elseif seg==3
		ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[mesh_collection.mesh.cell[:,i_mesh]][3]
	    a_2 = U[mesh_collection.mesh.cell[:,i_mesh]][1]

		basis_left = u0[ind_edge,3]
		basis_right = u0[ind_edge,1]
	end


	k_edge = 0.01

	system_matrix = [ (a_1^2 + k_edge)*speye(n)   (a_1*a_2)*speye(n) ; 
                        (a_1*a_2)*speye(n)   (a_2^2 + k_edge)*speye(n) ]

    rhs = ( [k_edge*basis_left + a_1*uOrigEdge ;
            k_edge*basis_right + a_2*uOrigEdge]
            - system_matrix[:,ind_con]*val_con )


    basis_lr = zeros(2*n)
    basis_lr[ind_con] = val_con
    basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]



    display([norm(uOrigEdge-a_1*basis_left - a_2*basis_right)  norm(uOrigEdge-a_1*basis_lr[1:n] - a_2*basis_lr[n+1:end])])

    plot(linspace(0,1,n), uOrigEdge, color="black", linewidth=2.5)

    plot(linspace(0,1,n), basis_lr[1:n], color="orange", linewidth=1.5)
    plot(linspace(0,1,n), basis_lr[n+1:end], color="orange", linewidth=1.5)

    plot(linspace(0,1,n), a_1*basis_lr[1:n] + a_2*basis_lr[n+1:end], color="green", linewidth=2.5)


	
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






