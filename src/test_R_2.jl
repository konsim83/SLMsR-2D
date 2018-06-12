# ---------------------------------------------------------------------------

import Parameter, Problem, FEM, Solver, PostProcess, Vis, Reconstruction

# ---------------------------------------------------------------------------


compute_low = true
compute_ref = true
compute_ms_reconstruction = true

post_process = true

# ---------------------------------------------------------------------------
# -------   Problem Parameters   -------
T_max = 1.0


# Multiscale diffusion, constant advection
# problem = Problem.Gaussian_R_1(T_max, 30)


# Constant low diffusion, solenoidal, traveling vortex
problem = Problem.Gaussian_R_2(T_max, 0.3 , 30)
# problem = Problem.Gaussian_R_2_conserv(T_max, 0.3 , 30)


# Multiscale diffusion, solenoidal, traveling vortex
# problem = Problem.Gaussian_R_3(T_max, 0.2 , 30)
# problem = Problem.Gaussian_R_3_conserv(T_max, 0.1 , 30)


# Constant low diffusion, divergent, traveling vortex
# problem = Problem.Gaussian_R_4(T_max, 0.2 , 30)
# problem = Problem.Gaussian_R_4_conserv(T_max, 0.1 , 30)


# Multiscale diffusion, divergent, traveling vortex
# problem = Problem.Gaussian_R_5(T_max, 0.2 , 30)
# problem = Problem.Gaussian_R_5_conserv(T_max, 0.1 , 30)
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# -------   Mesh parameters   -------
n_edge_per_seg = 4
n_refinement = 4


n_edge_per_seg_f = 15
max_are_cell_f = 0.005

n_edge_per_seg_f = 20
max_are_cell_f = 0.004

# n_edge_per_seg_f = 26
# max_are_cell_f = 0.001


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 3

n_order_FEM_f = 1
n_order_quad_f = n_order_quad


time_step_method = 1

dt = 1/100
n_steps_f = 5

# ---------------------------
# 1: non-conformal L2-reconstruction
# 2: conformal L2-reconstruction
# 3: conformal (smooth) H1-reconstruction
# 4: conformal (smooth) H1-reconstruction with soft/hard PoU constraint
reconstruction_method = 3

reconstruct_edge = false
# ---------------------------



k_edge = 0.01
k_int = [1.0 ; 1.0 ; 1.0] * 0.01
k_sum = 0
k = [k_int ; k_edge; k_sum]
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
if compute_low
        # -------   Build parameter structure   -------
        par_FEM_low = Parameter.Parameter_FEM(problem.T,
                                                dt,
                                                n_edge_per_seg,
                                                0,
                                                n_order_FEM,
                                                n_order_quad,
                                                time_step_method)

        # -------   Call the solver   -------
        @time solution_FEM_low, mesh_FEM_low = Solver.solve_FEM_periodic_square(par_FEM_low, problem)
end

if compute_ref
        # -------   Build parameter structure   -------
        par_FEM_ref = Parameter.Parameter_FEM(problem.T,
                                                dt,
                                                n_edge_per_seg,
                                                n_refinement,
                                                n_order_FEM,
                                                n_order_quad,
                                                time_step_method)

        # -------   Call the solver   -------
        @time solution_FEM_ref, mesh_FEM_ref = Solver.solve_FEM_periodic_square(par_FEM_ref, problem)
        PostProcess.evaluateGrad!(solution_FEM_ref, mesh_FEM_ref)
end


if compute_ms_reconstruction
        # -------   Build parameter structure   -------
        par_MsFEM_r = Parameter.Parameter_MsFEM(problem.T,
                                                dt,
                                                n_steps_f,
                                                n_edge_per_seg,
                                                n_refinement,
                                                n_edge_per_seg_f,
                                                max_are_cell_f,
                                                n_order_FEM_f,
                                                n_order_quad_f,
                                                time_step_method,
                                                reconstruction_method,
                                                reconstruct_edge,
                                                k)


        # -------   Call the solver   -------
        @time solution_ms_r, mesh_collection_r = Solver.solve_MsFEM_periodic_square_reconstruction(par_MsFEM_r, problem)
end
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
if post_process
        # ---------------------
        solution_FEM_low_mapped = PostProcess.map_solution(solution_FEM_low, mesh_FEM_low, mesh_FEM_ref)
        solution_ms_r_mapped = PostProcess.map_solution(solution_ms_r, mesh_collection_r, mesh_FEM_ref)
        # ---------------------
        

        # ---------------------
        err_low = PostProcess.error_L2(solution_FEM_ref, solution_FEM_low_mapped)[:]
        err_low_H1 = PostProcess.error_H1(solution_FEM_ref, solution_FEM_low_mapped)[:]
        Vis.writeError2file(err_low, T_max, "Error-L2-" * problem.file_name * "-FEM-low", "L2")
        Vis.writeError2file(err_low_H1, T_max, "Error-H1-" * problem.file_name * "-FEM-low", "H1")
        # ---------------------


        # ---------------------
        if reconstruct_edge
            err_ms_r_recon = PostProcess.error_L2(solution_FEM_ref, solution_ms_r_mapped)[:]
            err_ms_r_recon_H1 = PostProcess.error_H1(solution_FEM_ref, solution_ms_r_mapped)[:]

            Vis.writeError2file(err_ms_r_recon, T_max, "Error-L2-" * problem.file_name * "-MsFEM_r-reconstruct-mapped", "L2")
            Vis.writeError2file(err_ms_r_recon_H1, T_max, "Error-H1-" * problem.file_name * "-MsFEM_r-reconstruct-mapped", "H1")
        else
            err_ms_r = PostProcess.error_L2(solution_FEM_ref, solution_ms_r_mapped)[:]
            err_ms_r_H1 = PostProcess.error_H1(solution_FEM_ref, solution_ms_r_mapped)[:]

            Vis.writeError2file(err_ms_r, T_max, "Error-L2-" * problem.file_name * "-MsFEM_r-mapped", "L2")
            Vis.writeError2file(err_ms_r_H1, T_max, "Error-H1-" * problem.file_name * "-MsFEM_r-mapped", "H1")
        end
        # ---------------------


        # ---------------------
        Vis.writeSolution_all(solution_FEM_low, mesh_FEM_low, problem.file_name * "-FEM-low")

        Vis.writeSolution_all(solution_FEM_ref, mesh_FEM_ref, problem.file_name * "-FEM-ref")
        
        Vis.writeSolution_all(solution_FEM_low_mapped, mesh_FEM_ref, problem.file_name * "-FEM-low-mapped")
        # ---------------------
        
        # ---------------------
        if reconstruct_edge
            Vis.writeSolution_all(solution_ms_r_mapped, mesh_FEM_ref, problem.file_name * "-MsFEM_r-reconstruct-mapped")
            
            i,j = ind2sub(mesh_collection_r.mesh.cell,find(mesh_collection_r.mesh.cell.==35))
            for ind in j
               Vis.writeBasis_all_steps(solution_ms_r, mesh_collection_r, ind, "Basis---recon---" * problem.file_name)
            end
        else
            Vis.writeSolution_all(solution_ms_r_mapped, mesh_FEM_ref, problem.file_name * "-MsFEM_r-mapped")
            
            i,j = ind2sub(mesh_collection_r.mesh.cell,find(mesh_collection_r.mesh.cell.==35))
            for ind in j
               Vis.writeBasis_all_steps(solution_ms_r, mesh_collection_r, ind, "Basis---" * problem.file_name)
            end
        end
        # ---------------------
end
# ---------------------------------------------------------------------------