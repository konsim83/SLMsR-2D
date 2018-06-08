function evolve_edge!(uBasisOpt :: Array{Float64,2},,
						mesh_local  :: Mesh.TriangleMesh.TriMesh,						
						point_orig :: Array{Array{Float64,2},1},
						problem_local :: Problem.AbstractBasisProblem,
						dt_f :: Float64, n_steps_f :: Int,
 						seg :: Int,
 						T :: Float64)

	# Connectivity
	cell_edge = sort(mesh_local.segment[:,mesh_local.segment_marker.==seg],1)
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
	quad_1D = Quad.Quad_line(0, 0, 2)

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
            M_orig = FEM_1D.assemble_mass(mesh_old, quad_1D, ref_el_1d, dof_old, problem_local)

            R_next = FEM.assemble_reaction(mesh_next, 
            								quad_1d,
											ref_el_1d,
            								dof_next,                                             
                                            problem_local,
                                            T - (n_steps_back-(j-1))*dt_f)

            D_next = FEM_1D.assemble_diffusion(mesh_next, 
            								quad_1d,
											ref_el_1d,
            								dof_next,                                             
                                            problem_local,
                                            T - (n_steps_back-(j-1))*dt_f)

            # Zero forcing
	        f_orig = zeros(mesh_old.n_point,2)

			# ----------------------------------
			# -------   Mu_t + Ru = Du   -------
			system_matrix = M_orig - dt_f*(D_next-R_next)

			rhs = M_orig*uBasisOpt[ind_edge, ind_basis,j-1] + dt_f*f_orig - system_matrix[:,[1;n]]*eye(2)

	        uBasisOpt[ind_edge[2:end-1],ind_basis,j] = system_matrix[2:end-1,2:end-1] \ rhs[2:end-1]
	        uBasisOpt[ind_edge[1;n],ind_basis,j] = eye(2)
			# ----------------------------------
	    end # end for

	    return nothing
end