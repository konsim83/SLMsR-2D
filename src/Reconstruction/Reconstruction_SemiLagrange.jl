function SemiLagrange_L2_opt(solution :: FEM.Solution_MsFEM,
								timeStepper :: Array{TimeIntegrator.ImplEuler,1},
								mesh_collection :: Mesh.TriMesh_collection,
		                        dof_collection :: FEM.AbstractDofCollection,
		                        ref_el_f :: FEM.AbstractRefEl,
		                        quad_f :: Quad.AbstractQuad,
		                        par :: Parameter.Parameter_MsFEM,
		                        problem :: Problem.AbstractPhysicalProblem,
		                        problem_f :: Array{Problem.AbstractBasisProblem,1},
		                        T :: Float64, k_time :: Int)

	# Needs to be set up only once since velocity is the same everywhere
	velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

		return hcat(-Problem.velocity(problem, -t, x)...)
	end


	point_all = hcat([mesh.point for mesh in mesh_collection.mesh_f]...)
	@time point_orig_all = traceback(point_all, T, velocity, par)

	n_steps_back = length(point_orig_all)-1
    dt_f = par.dt/(n_steps_back)

    point_count = 0
	for i in 1:mesh_collection.mesh.n_cell
		mesh_local = mesh_collection.mesh_f[i]
		ind = (1:mesh_local.n_point) + point_count

		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Reconstruct the initial values
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		u_basis_tmp = zeros(mesh_local.n_point, 3, n_steps_back+1)
		u_basis_tmp[:,:,1] = Problem.u_init(problem_f[i], mesh_local.point)

		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		# ----------------------------------------------------
        # -------   Time evolution on distinct grids   -------
        for j=2:(n_steps_back+1)

        	p_old = point_orig_all[j-1][:,ind]
			p_next = point_orig_all[j][:,ind]

            # Create mesh of next time step
            mesh_old = Mesh.TriMesh(mesh_local, p_old, "Old mesh...")
            mesh_next = Mesh.TriMesh(mesh_local, p_next, "Next mesh...")

            # Create dof handler of next timestep
            dof_old = FEM.Dof_Pk(mesh_old, problem_f[i], 1)
            dof_next = FEM.Dof_Pk(mesh_next, problem_f[i], 1)

            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        # Time step (small) with implicit Euler
            M_orig = FEM.assemble_mass(mesh_old, dof_old, ref_el_f, quad_f, problem_f[i])
            M_next = FEM.assemble_mass(mesh_next, dof_next, ref_el_f, quad_f, problem_f[i])
            M = 0.5*(M_next + M_orig)

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
	        # Set the system matrices 
	        TimeIntegrator.updateSystem!(timeStepper[i].systemData, M, D_next-R_next, f_orig, dof_next, u_basis_tmp[:,:,j-1], par.dt)

	        # Make a single time step
	        TimeIntegrator.makeStep!(timeStepper[i], dof_next, view(u_basis_tmp, :, :,j), u_basis_tmp[:,:,j-1])
			# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
        end # end for
        # ----------------------------------------------------
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Pass the result of the evolution to the basis...
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		solution.phi_1[i][:,k_time+1] = u_basis_tmp[:,1,end]
        solution.phi_2[i][:,k_time+1] = u_basis_tmp[:,2,end]
        solution.phi_3[i][:,k_time+1] = u_basis_tmp[:,3,end]
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		point_count += mesh_local.n_point
	end # end for


	return point_orig_all
end


function traceback(point :: Array{Float64,2}, T :: Float64, 
					velocity :: Function, par :: Parameter.Parameter_MsFEM,)

	tspan = (-T, -T+par.dt)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob)

	return sol.u
end