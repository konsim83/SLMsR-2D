function SemiLagrange_L2_opt(solution :: FEM.Solution_MsFEM,
								mesh_collection :: Mesh.TriMesh_collection,
		                        dof_collection :: FEM.AbstractDofCollection,
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

	n_steps_back = length(point_orig)-1
    dt_f = par.dt/(n_steps_back)

    point_count = 0
	for i in 1:mesh_collection.mesh.n_cell
		mesh_local = mesh_collection.mesh_f[i]
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Reconstruct the initial values
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		# --------------------------------------------------------------
        # -------   Time evolution on distinct grids   -------
        ind = (1:mesh_local.n_point) + point_count
        for j=2:(n_steps_back+1)

        	p_old = point_orig_all[j-1][:,ind]
			p_next = = point_orig_all[j][:,ind]

            # Create mesh of next time step
            mesh_old = Mesh.TriangeMesh.TriMesh(mesh_local, p_old, "Old mesh...")
            mesh_next = Mesh.TriangeMesh.TriMesh(mesh_local, p_next, "Next mesh...")

            # Create dof handler of next timestep
            dof_old = FEM.Dof_Pk(mesh_old, problem_f[i], 1)
            dof_next = FEM.Dof_Pk(mesh_next, problem_f[i], 1)


        	# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Time step (small) with implicit Euler

            M_orig = FEM.assemble_mass(mesh_old, dof_old, solution.ref_el_f, solution.ref_el_f.quad, problem_f[i])
            M_next = FEM.assemble_mass(mesh_next, dof_next, solution.ref_el_f, solution.ref_el_f.quad, problem_f[i])
            M = 0.5*(M_next + M_orig)

            f_orig = zeros(size(M_next,1))








            R_next = FEM.assemble_global_R(dof_next, 
                                            mesh_next, 
                                            solution.ref_el_f, 
                                            solution.ref_el_f.quad, 
                                            problem_f[i], 
                                            T - (n_steps_back-(j-1))*dt_f)
            D_next = FEM.assemble_global_D(dof_next, 
                                            mesh_next, 
                                            solution.ref_el_f,
                                            solution.ref_el_f.quad, 
                                            problem_f[i], 
                                            T - (n_steps_back-(j-1))*dt_f)

            # -------   u_t + Ru = Du   -------
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
        end # end for
        # --------------------------------------------------------------
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		# Pass the result of the evolution to the basis...
		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




		# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		point_count += mesh_local.n_point
	end # end for


	return point_orig
end


function traceback(point :: Array{Float64,2}, T :: Float64, 
					velocity :: Function, par :: Parameter.Parameter_MsFEM,)

	tspan = (-T, -T+par.dt)
	prob = DifferentialEquations.ODEProblem(velocity, point, tspan)
	sol = DifferentialEquations.solve(prob)

	return sol.u
end