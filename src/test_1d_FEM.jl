import Mesh, Problem, FEM, FEM_1D, Quad, TimeIntegrator, Plot

import DifferentialEquations

using PyPlot
PyPlot.pyimport("mpl_toolkits.mplot3d.axes3d")


smooth = false
# smooth = true


function traceback(point :: Array{Float64,2},
					T :: Float64, 
					velocity :: Function)
	dt = 1/10
	n_steps_back = 5

	tspan = (-T, -T+dt)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob, dtmax=dt/n_steps_back)

	return sol.u
end

# ----------------------------------------------------------------
m_coarse = Mesh.mesh_unit_square(4)
m_simplex = Mesh.mesh_unit_simplex(12, 0.005)

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

T = 1.0
i_mesh = 5 # mesh index

# problem = Problem.Gaussian_R_4(1.0, 0.2 , 30)
problem = Problem.Gaussian_R_4_conserv(1.0, 0.2 , 30)


# ----------------------------------------------------------------
velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

	return hcat(-Problem.velocity(problem, -t, x)...)
end


close()
close()
close()
close()
close()

mesh = mesh_collection.mesh_f[i_mesh]

tri = Geometry.Triangle(mesh_collection.mesh.point[:, mesh_collection.mesh.cell[:,i_mesh]])
problem_f = Problem.BasisFun(problem, tri)


	
point_orig = traceback(mesh.point, T, velocity)
pOrig = point_orig[end]
u0 = Problem.u_init(problem_f, point_orig[1])
uOrig = Problem.u_init(problem, pOrig)

n_steps_back = length(point_orig)-1
dt_f = (1/10)/n_steps_back

uBasisOpt = zeros(size(u0,1), size(u0,2), n_steps_back+1)

# -----------------------------------------
ind_corner = [sort(circshift(mesh.segment[:,mesh.segment_marker.==1],1)[:])[1] ; 
                sort(circshift(mesh.segment[:,mesh.segment_marker.==1],1)[:])[end] ; 
                sort(circshift(mesh.segment[:,mesh.segment_marker.==2],1)[:])[end] ]
pOrigCorner = pOrig[:,ind_corner]
U = Problem.u_init(problem, pOrigCorner)
# -----------------------------------------
	
# ----------------------------------------------------------------
for seg in 1:3

	cell = circshift(mesh.segment[:,mesh.segment_marker.==seg],1)
	ind_edge = unique(cell)
	n = length(ind_edge)	

	pOrigEdge = pOrig[:,ind_edge]
	uOrigEdge = uOrig[ind_edge]
	
	if seg==1
		ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[1]
	    a_2 = U[2]

		basis_left = u0[ind_edge,1]
		basis_right = u0[ind_edge,2]
	elseif seg==2
		ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[2]
	    a_2 = U[3]

		basis_left = u0[ind_edge,2]
		basis_right = u0[ind_edge,3]
	elseif seg==3
		ind_con = [1 ; n ; 1+n ; 2*n]
	    val_con = [1. ; 0. ; 0. ; 1.]
	    ind_uncon = setdiff(1:(2*n), ind_con)

	    a_1 = U[3]
	    a_2 = U[1]

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

    if seg==1
    	uBasisOpt[ind_edge,1,1] = basis_lr[1:n]
    	uBasisOpt[ind_edge,2,1] = basis_lr[n+1:end]
	elseif seg==2
		uBasisOpt[ind_edge,2,1] = basis_lr[1:n]
		uBasisOpt[ind_edge,3,1] = basis_lr[n+1:end]
	elseif seg==3
		uBasisOpt[ind_edge,3,1] = basis_lr[1:n]
		uBasisOpt[ind_edge,1,1] = basis_lr[n+1:end]
	end

	f1 = figure(1, figsize = (18,6))
	subplot(1,3,seg)
    plot(linspace(0,1,n), uOrigEdge, color="black", linewidth=2.5)
    plot(linspace(0,1,n), basis_lr[1:n], color="orange", linewidth=1.5)
    plot(linspace(0,1,n), basis_lr[n+1:end], color="orange", linewidth=1.5)
    plot(linspace(0,1,n), a_1*basis_lr[1:n] + a_2*basis_lr[n+1:end], color="green", linewidth=2.5)
end # end for seg



R = []
D = []


for seg in 1:3

	cell_edge = sort(mesh.segment[:,mesh.segment_marker.==seg],1)
	if seg==3
		cell_edge[[1;2],end] = cell_edge[[2;1],end]
	end

	if seg==1
		ind_basis = [1;2]
	elseif seg==2
		ind_basis = [2;3]
	elseif seg==3
		ind_basis = [3;1]
	end

    # -----------------------     
    # We need this array to plug the solution of
    # the evolution into the uBasisOpt vector.
	ind_edge = unique(cell_edge)
	n = length(ind_edge)
	# -----------------------

	ref_el_1d = FEM_1D.RefEl_Lagrange_1()
	quad_1d = Quad.Quad_line(0, 0, 2)


	for j=2:(n_steps_back+1)
    	p_old = point_orig[end-(j-2)][:,ind_edge]
		p_next = point_orig[end-(j-1)][:,ind_edge]

        # Create mesh of next time step
        mesh_old = FEM_1D.Mesh_1D(cell_edge, p_old)
        mesh_next = FEM_1D.Mesh_1D(cell_edge, p_next)

        # Create dof handler of next timestep
        dof_old = FEM_1D.Dof_1D(mesh_old)
        dof_next = FEM_1D.Dof_1D(mesh_next)

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Time step (small) with implicit Euler
        M_orig = FEM_1D.assemble_mass(mesh_old, quad_1d, ref_el_1d, dof_old, problem_f)

        if problem_f.conservative
            R_next = FEM_1D.assemble_reaction(mesh_next, 
            								quad_1d,
											ref_el_1d,
            								dof_next,                                             
                                            problem_f,
                                            T - (n_steps_back-(j-1))*dt_f)
        else
        	R_next = spzeros(size(M_orig,1), size(M_orig,2))
        end

        D_next = FEM_1D.assemble_diffusion(mesh_next, 
            								quad_1d,
											ref_el_1d,
            								dof_next,                                             
                                            problem_f,
                                            T - (n_steps_back-(j-1))*dt_f)

        R = R_next
        D = D_next
        # display((R_next))
        # display(full(D_next))

        # Zero forcing
        f_orig = zeros(mesh_old.n_point,2)

		# ----------------------------------
		# -------   Mu_t + Ru = Du   -------
		system_matrix = M_orig - dt_f*(D_next-R_next)

		rhs = M_orig*uBasisOpt[ind_edge, ind_basis,j-1] + dt_f*f_orig - system_matrix[:,[1;n]]*eye(2)

        uBasisOpt[ind_edge[2:end-1], ind_basis, j] = system_matrix[2:end-1,2:end-1] \ rhs[2:end-1,:]
        uBasisOpt[ind_edge[[1;n]], ind_basis, j] = eye(2)
		# ----------------------------------
    end # end for


end # end for evolve seg















#=
# -----------------------------
n_dof = mesh.n_point
nOrig = length(uOrig)

u1 = U[1]
u2 = U[2]
u3 = U[3]

uOpt = vec(uBasisOpt)

ind_con = [find(mesh.point_marker.!=0) ; 
			find(mesh.point_marker.!=0)+n_dof ;
			find(mesh.point_marker.!=0)+2*n_dof]			
constr_val = uOpt[ind_con]
ind_uncon = setdiff(1:(3*nOrig), ind_con)

nOrig = length(uOrig)

if smooth
	ref_el = FEM.RefEl_Pk(1)
    quad = Quad.Quad_simplex(3)
    dof = FEM.Dof_Pk(mesh, problem, 1)

    k_smooth = 0.01

	laplace_matrix = FEM.assemble_Laplace(mesh, 
                                            dof,
                                            ref_el, 
                                            quad)
    # A = laplace_matrix' * laplace_matrix
    A = laplace_matrix

    system_matrix = [u1^2*speye(nOrig)+k_smooth*A    u1*u2*speye(nOrig)               u1*u3*speye(nOrig) ; 
                     u2*u1*speye(nOrig)              u2^2*speye(nOrig)+k_smooth*A     u2*u3*speye(nOrig) ; 
                     u3*u1*speye(nOrig)              u3*u2*speye(nOrig)               u3^2*speye(nOrig)+k_smooth*A ]

    rhs = [ u1*uOrig ;
            u2*uOrig ;
            u3*uOrig ] - system_matrix[:,ind_con]*constr_val
else
	k = [1 ; 1 ; 1] * 0.01

	system_matrix = [(u1^2+k[1])*speye(nOrig)    u1*u2*speye(nOrig)          u1*u3*speye(nOrig) ; 
                    u2*u1*speye(nOrig)           (u2^2+k[2])*speye(nOrig)    u2*u3*speye(nOrig) ; 
                    u3*u1*speye(nOrig)           u3*u2*speye(nOrig)          (u3^2+k[3])*speye(nOrig)]

	rhs = [k[1]*u0[:,1] + u1*uOrig; 
	        k[2]*u0[:,2] + u2*uOrig ; 
	        k[3]*u0[:,3] + u3*uOrig]  - system_matrix[:,ind_con]*constr_val
end

uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]

uOpt = reshape(uOpt,:,3)

println("\n---------- Surface approximation ----------")
display([norm(uOrig-(u1*u0[:,1]+u2*u0[:,2]+u3*u0[:,3]))   norm(uOrig-(u1*uOpt[:,1]+u2*uOpt[:,2]+u3*uOpt[:,3]))])

f2 = figure(2, figsize = (18,6))
ax1 = subplot(131, projection="3d")
ax1[:plot_trisurf](pOrig[1,:], pOrig[2,:], uOrig, triangles=mesh.cell'-1)
ax2 = subplot(132, projection="3d")
ax2[:plot_trisurf](pOrig[1,:], pOrig[2,:], u1*u0[:,1]+u2*u0[:,2]+u3*u0[:,3], triangles=mesh.cell'-1)
ax3 = subplot(133, projection="3d")
ax3[:plot_trisurf](pOrig[1,:], pOrig[2,:], u1*uOpt[:,1]+u2*uOpt[:,2]+u3*uOpt[:,3], triangles=mesh.cell'-1)

f3 = figure(3, figsize = (18,6))
subplot(131, projection="3d")
plot_trisurf(pOrig[1,:], pOrig[2,:], uOpt[:,1], triangles=mesh.cell'-1)
subplot(132, projection="3d")
plot_trisurf(pOrig[1,:], pOrig[2,:], uOpt[:,2], triangles=mesh.cell'-1)
subplot(133, projection="3d")
plot_trisurf(pOrig[1,:], pOrig[2,:], uOpt[:,3], triangles=mesh.cell'-1)
=#