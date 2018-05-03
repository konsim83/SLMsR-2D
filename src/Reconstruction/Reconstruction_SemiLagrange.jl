function reconstruct_L2(solution :: FEM.Solution_MsFEM,
							mesh_collection :: Mesh.TriMesh_collection,
							par :: Parameter.Parameter_MsFEM,
							problem_f :: Problem.AbstractBasisProblem,
							point_orig :: Array{Float64,2},
							u_orig :: Array{Float64,1},
							k_time :: Int,
							ind_cell :: Int)

    # ------------------------------------------------
    # Evaluate the solution at the traced back points
    # u = PostProcess.evaluate(solution, mesh_collection, point_orig, k_time-1)
    u = u_orig
    # ------------------------------------------------

	u0 = Problem.u_init(problem_f, mesh_collection.mesh_f[ind_cell].point)

    k = par.k

    # Contraints for nodal values of basis
    n_dof = mesh_collection.mesh_f[ind_cell].n_point
    ind_con = [1 ; 2 ; 3 ; n_dof+1 ; n_dof+2 ; n_dof+3 ; 2*n_dof+1 ; 2*n_dof+2 ; 2*n_dof+3]
    constr_val = vec(eye(3))
    ind_uncon = setdiff(1:(3*n_dof), ind_con)

    U = solution.u[mesh_collection.mesh.cell[:,ind_cell],k_time-1]
    u1 = U[1]
    u2 = U[2]
    u3 = U[3]

    n = length(u)
    system_matrix = [(u1^2+k[1])*speye(n) u1*u2*speye(n) u1*u3*speye(n) ; 
    					u2*u1*speye(n) (u2^2+k[2])*speye(n) u2*u3*speye(n) ; 
    					u3*u1*speye(n) u3*u2*speye(n) (u3^2+k[3])*speye(n)]
	rhs = [k[1]*u0[:,1] + u1*u; 
			k[2]*u0[:,2] + u2*u ; 
			k[3]*u0[:,3] + u3*u]  - system_matrix[:,ind_con]*constr_val


	uOpt = zeros(3*n_dof)
	uOpt[ind_con] = constr_val
	uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]

	return reshape(uOpt,:,3)
end


function reconstruct_L2(solution :: FEM.Solution_MsFEM,
							mesh_collection :: Mesh.TriMesh_collection,
							par :: Parameter.Parameter_MsFEM,
							problem :: Problem.AbstractPhysicalProblem,
							problem_f :: Problem.AbstractBasisProblem,
							point :: Array{Float64,2},
							ind_cell :: Int)
    
    # ------------------------------------------------
    # Evaluate the solution at the traced back points
    u = Problem.u_init(problem, point)
    # ------------------------------------------------

	u0 = Problem.u_init(problem_f, mesh_collection.mesh_f[ind_cell].point)

    k = par.k

    # Contraints for nodal values of basis
    n_dof = mesh_collection.mesh_f[ind_cell].n_point
    ind_con = [1 ; 2 ; 3 ; n_dof+1 ; n_dof+2 ; n_dof+3 ; 2*n_dof+1 ; 2*n_dof+2 ; 2*n_dof+3]
    constr_val = vec(eye(3))
    ind_uncon = setdiff(1:(3*n_dof), ind_con)

    # These are the weights of the basis functions
    U = solution.u[mesh_collection.mesh.cell[:,ind_cell],1]
    u1 = U[1]
    u2 = U[2]
    u3 = U[3]

    n = length(u)
    system_matrix = [(u1^2+k[1])*speye(n) u1*u2*speye(n) u1*u3*speye(n) ; 
    					u2*u1*speye(n) (u2^2+k[2])*speye(n) u2*u3*speye(n) ; 
    					u3*u1*speye(n) u3*u2*speye(n) (u3^2+k[3])*speye(n)]
	rhs = [k[1]*u0[:,1] + u1*u; 
			k[2]*u0[:,2] + u2*u ; 
			k[3]*u0[:,3] + u3*u] - system_matrix[:,ind_con]*constr_val


	uOpt = zeros(3*n_dof)
	uOpt[ind_con] = constr_val
	uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]

	return reshape(uOpt,:,3)
end


function SemiLagrange_L2_opt!(solution :: FEM.Solution_MsFEM,
								timeStepper :: Array{TimeIntegrator.ImplEuler,1},
								mesh_collection :: Mesh.TriMesh_collection,
		                        dof_collection :: FEM.AbstractDofCollection,
		                        ref_el_f :: FEM.AbstractRefEl,
		                        quad_f :: Quad.AbstractQuad,
		                        par :: Parameter.Parameter_MsFEM,
		                        problem :: Problem.AbstractPhysicalProblem,
		                        problem_f :: Array{Problem.AbstractBasisProblem,1},
		                        T :: Float64, k_time :: Int)

	# Reconstruct the first basis as well but only at the beginning, i.e.,
	# k_time==2.
	if k_time==2
		for i in 1:mesh_collection.mesh.n_cell
			u_basis_tmp = zeros(mesh_collection.mesh_f[i].n_point, 3)
			u_basis_tmp[:,:] = reconstruct_L2(solution,
													mesh_collection,
													par,
													problem,
													problem_f[i],
													mesh_collection.mesh_f[i].point,
													i)
			solution.phi_1[i][:,k_time-1] = u_basis_tmp[:,1]
		    solution.phi_2[i][:,k_time-1] = u_basis_tmp[:,2]
		    solution.phi_3[i][:,k_time-1] = u_basis_tmp[:,3]
		end
	end

	# Needs to be set up only once since velocity is the same everywhere
	velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

		return hcat(-Problem.velocity(problem, -t, x)...)
	end


	point_all = hcat([mesh.point for mesh in mesh_collection.mesh_f]...)

	display("---------------------------------------")
	display("Reconstruction at time step $(k_time) / $(par.n_steps+1):")
	point_orig_all = traceback(point_all, T, velocity, par)

	if k_time==2
		@time u_orig_all = Problem.u_init(problem, point_orig_all[end])
	else
		@time u_orig_all = PostProcess.evaluate(solution, mesh_collection, point_orig_all[end], k_time-1)
	end

	n_steps_back = length(point_orig_all)-1
    dt_f = par.dt/n_steps_back

    point_count = 0
	@time for i in 1:mesh_collection.mesh.n_cell
		mesh_local = mesh_collection.mesh_f[i]
		ind = (1:mesh_local.n_point) + point_count

		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Reconstruct the initial values
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		u_basis_tmp = zeros(mesh_local.n_point, 3, n_steps_back+1)
		# u_basis_tmp[:,:,1] = Problem.u_init(problem_f[i], mesh_local.point)
		u_basis_tmp[:,:,1] = reconstruct_L2(solution,
												mesh_collection,
												par,
												problem_f[i],
												point_orig_all[end][:,ind],
												u_orig_all[ind],
												k_time,
												i)

		# Contraints for nodal values of basis
	    n_dof = mesh_local.n_point
	    ind_con = [1 ; 2 ; 3]
	    constr_val = eye(3)
	    ind_uncon = setdiff(1:n_dof, ind_con)

		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Evolve the reconstructed boundaries if for reaction
		# problems
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if problem.conservative
			# Solve transient 2D problems
			for i=1:3
				# Indices of points on certain boundary edge (indices in terms
				# of original local mesh)
				# cell_2d = mesh_local.segment[:,mesh_local.segment_marker.==i] # 2-by-n matrix


	        
	        	u_basis_tmp[4:end,:,j] = 

			end
		end

		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		# ----------------------------------------------------
        # -------   Time evolution on distinct grids   -------
        for j=2:(n_steps_back+1)
        	p_old = point_orig_all[end-(j-2)][:,ind]
			p_next = point_orig_all[end-(j-1)][:,ind]

            # Create mesh of next time step
            mesh_old = Mesh.TriMesh(mesh_local, p_old, "Old mesh...")
            mesh_next = Mesh.TriMesh(mesh_local, p_next, "Next mesh...")

            # Create dof handler of next timestep
            dof_old = FEM.Dof_Pk(mesh_old, problem_f[i], 1)
            dof_next = FEM.Dof_Pk(mesh_next, problem_f[i], 1)

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        # Time step (small) with implicit Euler
            M_orig = FEM.assemble_mass(mesh_old, dof_old, ref_el_f, quad_f, problem_f[i])
            # M_next = FEM.assemble_mass(mesh_next, dof_next, ref_el_f, quad_f, problem_f[i])
            # M = 0.5*(M_next + M_orig)

            R_next = FEM.assemble_reaction(mesh_next, 
            								dof_next,
            								ref_el_f, 
                                            quad_f, 
                                            problem_f[i],
                                            T - (n_steps_back-(j-1))*dt_f)

            D_next = FEM.assemble_diffusion(mesh_next, 
            								dof_next,
            								ref_el_f, 
                                            quad_f, 
                                            problem_f[i],
                                            T - (n_steps_back-(j-1))*dt_f)

            # Zero forcing
	        f_orig = zeros(mesh_local.n_point,3)

			# -------   u_t + Ru = Du   -------
	        if problem.conservative
	        	# system_matrix = M_orig-dt_f*(D_next-R_next)
	        	# rhs = M_orig*u_basis_tmp[:,:,j-1] - system_matrix[:,ind_con]*constr_val
	        	# u_basis_tmp[ind_con,:,j] = constr_val
	        	# u_basis_tmp[ind_uncon,:,j] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon,:]

	        	# Set the system matrices, BV from new state
		        TimeIntegrator.updateSystem!(timeStepper[i].systemData, M_orig, D_next-R_next, f_orig, dof_next, u_basis_tmp[:,:,j], dt_f)

		        # Make a single time step
		        TimeIntegrator.makeStep!(timeStepper[i], dof_next, view(u_basis_tmp, :, :,j), u_basis_tmp[:,:,j-1])
	        else
		        # Set the system matrices 
		        TimeIntegrator.updateSystem!(timeStepper[i].systemData, M_orig, D_next-R_next, f_orig, dof_next, u_basis_tmp[:,:,j-1], dt_f)

		        # Make a single time step
		        TimeIntegrator.makeStep!(timeStepper[i], dof_next, view(u_basis_tmp, :, :,j), u_basis_tmp[:,:,j-1])
		    end
			# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
        end # end for
        # ----------------------------------------------------
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Pass the result of the evolution to the basis...
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		solution.phi_1[i][:,k_time] = u_basis_tmp[:,1,end]
        solution.phi_2[i][:,k_time] = u_basis_tmp[:,2,end]
        solution.phi_3[i][:,k_time] = u_basis_tmp[:,3,end]
        
        # Compute time derivative of basis functions via finite differencing
        solution.phi_1_t[i][:,k_time] = FiniteDiff.backward_single(solution.phi_1[i][:,(k_time-1):k_time], par.dt)
        solution.phi_2_t[i][:,k_time] = FiniteDiff.backward_single(solution.phi_2[i][:,(k_time-1):k_time], par.dt)
        solution.phi_3_t[i][:,k_time] = FiniteDiff.backward_single(solution.phi_3[i][:,(k_time-1):k_time], par.dt)
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		point_count += mesh_local.n_point
	end # end for
	display("---------------------------------------")

	return nothing
end


function traceback(point :: Array{Float64,2}, T :: Float64, 
					velocity :: Function, par :: Parameter.Parameter_MsFEM)

	tspan = (-T, -T+par.dt)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob, dtmax=par.dt/par.n_steps_f)

	return sol.u
end


function traceback(point :: Array{Float64,2},
					T0 :: Float64, 
					Tend :: Float64, 
					velocity :: Function)

	tspan = (T0, Tend)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob)

	return sol.u
end