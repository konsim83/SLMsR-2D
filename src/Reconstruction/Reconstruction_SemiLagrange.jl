function SemiLagrange_opt!(solution :: FEM.Solution_MsFEM,
							timeStepper :: Array{TimeIntegrator.ImplEuler,1},
							mesh_collection :: Mesh.TriMesh_collection,
	                        dof_collection :: FEM.AbstractDofCollection,
	                        ref_el_f :: FEM.AbstractRefEl,
	                        quad_f :: Quad.AbstractQuad,
	                        par :: Parameter.Parameter_MsFEM,
	                        problem :: Problem.AbstractPhysicalProblem,
	                        problem_f :: Array{Problem.AbstractBasisProblem,1},
	                        T :: Float64,
	                        k_time :: Int)


	if k_time==2
		println("\n---------------------------------------")
		println("Reconstruction of initial basis:")
		@time for i_mesh_local in 1:mesh_collection.mesh.n_cell
			u_basis_tmp = zeros(mesh_collection.mesh_f[i_mesh_local].n_point, 3)
			if par.reconstruction_method==0
				u_basis_tmp[:,:] = reconstruct_L2_noConstraint(solution,
																mesh_collection,
																par,
																problem,
																problem_f[i_mesh_local],
																mesh_collection.mesh_f[i_mesh_local].point,
																i_mesh_local)
			elseif par.reconstruction_method==1
				u_basis_tmp[:,:] = reconstruct_L2_nonconformal(solution,
																mesh_collection,
																par,
																problem,
																problem_f[i_mesh_local],
																mesh_collection.mesh_f[i_mesh_local].point,
																i_mesh_local)
			elseif par.reconstruction_method==2
				u_basis_tmp[:,:] = reconstruct_L2_conformal(solution,
															mesh_collection,
															par,
															problem,
															problem_f[i_mesh_local],
															mesh_collection.mesh_f[i_mesh_local].point,
															i_mesh_local)
			elseif par.reconstruction_method==3
				u_basis_tmp[:,:] = reconstruct_H1_conformal(solution,
																mesh_collection,
																dof_collection,
																ref_el_f,
				                        						quad_f,
																par,
																problem,
																problem_f[i_mesh_local],
																mesh_collection.mesh_f[i_mesh_local].point,
																i_mesh_local)
			elseif par.reconstruction_method==4
				u_basis_tmp[:,:] = reconstruct_H1_conformal_sumEqualOne(solution,
																mesh_collection,
																dof_collection,
																ref_el_f,
				                        						quad_f,
																par,
																problem,
																problem_f[i_mesh_local],
																mesh_collection.mesh_f[i_mesh_local].point,
																i_mesh_local)
			end
			solution.phi_1[i_mesh_local][:,1] = u_basis_tmp[:,1]
		    solution.phi_2[i_mesh_local][:,1] = u_basis_tmp[:,2]
		    solution.phi_3[i_mesh_local][:,1] = u_basis_tmp[:,3]
		end
		println("---------------------------------------\n")
	end


	# ------------------------------------------------------------------
	# Needs to be set up only once since velocity is the same everywhere
	velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

		return hcat(-Problem.velocity(problem, -t, x)...)
	end
	# ------------------------------------------------------------------

	
	# ------------------------------------------------------------------
	if problem.conservative # || problem.reconstructEdge
		println("\n---------------------------------------")
		println("Reconstruction at time step $(k_time) / $(par.n_steps+1):")
		println("Edge evolution --- yes")
	else
		println("\n---------------------------------------")
		println("Reconstruction at time step $(k_time) / $(par.n_steps+1):")
		println("Edge evolution --- no")
	end
	@time for i_mesh_local in 1:mesh_collection.mesh.n_cell
		mesh_local = mesh_collection.mesh_f[i_mesh_local]
		problem_local = problem_f[i_mesh_local]

		# Trace back points
		point_orig = traceback(mesh_local.point, T, velocity, par)

		n_steps_back = length(point_orig)-1
    	dt_f = par.dt/n_steps_back

		# Reconstruct basis at previous time step
		if k_time==2
			uOrig = Problem.u_init(problem, point_orig[end])
		else
			uOrig = PostProcess.evaluate(solution, mesh_collection, point_orig[end], k_time-1)
		end

    	# initialize basis
    	u_basis_tmp = zeros(mesh_local.n_point, 3, n_steps_back+1)

    	# Solve inverse problem for basis
    	if par.reconstruction_method==0
    		u_basis_tmp[:,:,1] = reconstruct_L2_noConstraint(solution,
																mesh_collection,
																par,
																problem_local,
																uOrig,
																k_time,
																i_mesh_local)
    	elseif par.reconstruction_method==1
    		u_basis_tmp[:,:,1] = reconstruct_L2_nonconformal(solution,
																mesh_collection,
																par,
																problem_local,
																uOrig,
																k_time,
																i_mesh_local)
    	elseif par.reconstruction_method==2
    		u_basis_tmp[:,:,1] = reconstruct_L2_conformal(solution,
															mesh_collection,
															par,
															problem_local,
															uOrig,
															k_time,
															i_mesh_local)
		elseif par.reconstruction_method==3
	    	u_basis_tmp[:,:,1] = reconstruct_H1_conformal(solution,
															mesh_collection,
															dof_collection,
															ref_el_f,
			                        						quad_f,
															par,
															problem_local,
															uOrig,
															k_time,
															i_mesh_local)
    	elseif par.reconstruction_method==4
    		u_basis_tmp[:,:,1] = reconstruct_H1_conformal_sumEqualOne(solution,
															mesh_collection,
															dof_collection,
															ref_el_f,
			                        						quad_f,
															par,
															problem_local,
															uOrig,
															k_time,
															i_mesh_local)
	    end


    	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	# Evolve boundary values
		if problem.conservative  || par.reconstructEdge
			@time for segment in 1:3
				evolve_edge!(u_basis_tmp, mesh_local, point_orig, problem_local, dt_f, n_steps_back, segment, T)
			end
		end
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # -------   Time evolution on distinct grids   -------
        for j=2:(n_steps_back+1)
        	p_old = point_orig[end-(j-2)]
			p_next = point_orig[end-(j-1)]

            # Create mesh of next time step
            mesh_old = Mesh.TriMesh(mesh_local, p_old, "Old mesh...")
            mesh_next = Mesh.TriMesh(mesh_local, p_next, "Next mesh...")

            # Create dof handler of next timestep
            dof_old = FEM.Dof_Pk(mesh_old, problem_local, 1)
            dof_next = FEM.Dof_Pk(mesh_next, problem_local, 1)

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        # Time step (small) with implicit Euler
            M_orig = FEM.assemble_mass(mesh_old, dof_old, ref_el_f, quad_f, problem_local)
            # M_next = FEM.assemble_mass(mesh_next, dof_next, ref_el_f, quad_f, problem_local)
            # M = 0.5*(M_next + M_orig)

            R_next = FEM.assemble_reaction(mesh_next, 
            								dof_next,
            								ref_el_f, 
                                            quad_f, 
                                            problem_local,
                                            T - (n_steps_back-(j-1))*dt_f)

            D_next = FEM.assemble_diffusion(mesh_next, 
            								dof_next,
            								ref_el_f, 
                                            quad_f, 
                                            problem_local,
                                            T - (n_steps_back-(j-1))*dt_f)

            # Zero forcing
	        f_orig = zeros(mesh_local.n_point,3)

			# -------   Mu_t + Ru = Du   -------
			# Set the system matrices 
	      #   if problem.conservative # || problem.reconstructEdge
	      #   	# error("Not implemented")
	      #   	ind_con = [1 ; 2 ; 3]
			    # constr_val = eye(3)
			    # ind_uncon = setdiff(1:mesh_local.n_point, ind_con)

	      #   	system = M_orig-dt_f*(D_next-R_next)
	      #   	rhs = M_orig[ind_uncon,:]*u_basis_tmp[:,:,j-1] - system[ind_uncon,ind_con]*constr_val
	      #   	u_basis_tmp[ind_uncon,:,j] = system[ind_uncon,ind_uncon]\ rhs
	      #   	u_basis_tmp[ind_con,:,j] = constr_val
	      #   else
		        # # Set the system matrices 
		        TimeIntegrator.updateSystem!(timeStepper[i_mesh_local].systemData, M_orig, D_next-R_next, f_orig, 
		        								dof_next, u_basis_tmp[:,:,j-1], u_basis_tmp[:,:,j-1], dt_f)
		        # Make a single time step
		        TimeIntegrator.makeStep!(timeStepper[i_mesh_local], dof_next, view(u_basis_tmp, :, :,j), u_basis_tmp[:,:,j-1])
	    	# end
	    end # end for
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		# Pass the result of the evolution to the basis
		solution.phi_1[i_mesh_local][:,k_time] = u_basis_tmp[:,1,end]
        solution.phi_2[i_mesh_local][:,k_time] = u_basis_tmp[:,2,end]
        solution.phi_3[i_mesh_local][:,k_time] = u_basis_tmp[:,3,end]
        
        # Compute time derivative of basis functions via finite differencing
        solution.phi_1_t[i_mesh_local][:,k_time] = FiniteDiff.backward_single(solution.phi_1[i_mesh_local][:,(k_time-1):k_time], par.dt)
        solution.phi_2_t[i_mesh_local][:,k_time] = FiniteDiff.backward_single(solution.phi_2[i_mesh_local][:,(k_time-1):k_time], par.dt)
        solution.phi_3_t[i_mesh_local][:,k_time] = FiniteDiff.backward_single(solution.phi_3[i_mesh_local][:,(k_time-1):k_time], par.dt)
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
    end # end for
    println("---------------------------------------\n\n")
    # ----------------------------------------------------

	return nothing
end


function traceback(point :: Array{Float64,2}, T :: Float64, 
					velocity :: Function, par :: Parameter.Parameter_MsFEM)

	tspan = (-T, -T+par.dt)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob, dtmax=par.dt/par.n_steps_f)

	return sol.u
end