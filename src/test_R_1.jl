# ---------------------------------------------------------------------------

import Parameter, Problem, FEM, Solver, PostProcess, Vis, Reconstruction

# ---------------------------------------------------------------------------


computeSTD = true
computeREF = true
computeMSR = true

post_process = true

# ---------------------------------------------------------------------------
# -------   Problem Parameters   -------
T_max = 1.0


# Multiscale diffusion, constant advection
problem = Problem.Gaussian_R_1(T_max, 30)


# Constant low diffusion, solenoidal, traveling vortex
# problem = Problem.Gaussian_R_2(T_max, 0.3 , 30)
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
n_refinement = 5


n_edge_per_seg_f = 15
max_area_cell_f = 0.005

n_edge_per_seg_f = 20
max_area_cell_f = 0.004

n_edge_per_seg_f = 26
max_area_cell_f = 0.001


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
if computeSTD
        # -------   Build parameter structure   -------
        parSTD = Parameter.Parameter_FEM(problem.T,
                                                dt,
                                                n_edge_per_seg,
                                                0, # no refinements
                                                n_order_FEM,
                                                n_order_quad,
                                                time_step_method)

        # -------   Call the solver   -------
        @time solSTD, meshSTD = Solver.solve_FEM_periodic_square(parSTD, problem)
end

if computeREF
        # -------   Build parameter structure   -------
        parREF = Parameter.Parameter_FEM(problem.T,
                                                dt,
                                                n_edge_per_seg,
                                                n_refinement,
                                                n_order_FEM,
                                                n_order_quad,
                                                time_step_method)

        # -------   Call the solver   -------
        @time solREF, meshREF = Solver.solve_FEM_periodic_square(parREF, problem)
        PostProcess.evaluateGrad!(solREF, meshREF)
end


if computeMSR
        # -------   Build parameter structure   -------
        parMSR = Parameter.Parameter_MsFEM(problem.T,
                                                dt,
                                                n_steps_f,
                                                n_edge_per_seg,
                                                n_refinement,
                                                n_edge_per_seg_f,
                                                max_area_cell_f,
                                                n_order_FEM_f,
                                                n_order_quad_f,
                                                time_step_method,
                                                reconstruction_method,
                                                reconstruct_edge,
                                                k)


        # -------   Call the solver   -------
        @time solMSR, meshMSR = Solver.solve_MsFEM_periodic_square_reconstruction(parMSR, problem)
end
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
if post_process
        # ---------------------
        solSTD_mapped = PostProcess.map_solution(solSTD, meshSTD, meshREF)
        solMSR_mapped = PostProcess.map_solution(solMSR, meshMSR, meshREF)
        # ---------------------
        

        # ---------------------
        # Compute and write standard error
        errSTD = PostProcess.error_L2(solREF, solSTD_mapped)[:]
        errSTD_H1 = PostProcess.error_H1(solREF, solSTD_mapped)[:]
        Vis.writeError2file(errSTD, T_max, "Error-L2-" * problem.file_name * "-STD", "L2")
        Vis.writeError2file(errSTD_H1, T_max, "Error-H1-" * problem.file_name * "-STD", "H1")
        # ---------------------


        # ---------------------
        # Write standard solutions
        Vis.writeSolution_all(solREF, meshREF, problem.file_name * "-REF")        
        Vis.writeSolution_all(solSTD_mapped, meshREF, problem.file_name * "-STD")
        # ---------------------


        # ---------------------
        # Compute and write multiscale errors and solutions
        if reconstruct_edge
            errMSR_recon = PostProcess.error_L2(solREF, solMSR_mapped)[:]
            errMSR_recon_H1 = PostProcess.error_H1(solREF, solMSR_mapped)[:]
            Vis.writeError2file(errMSR_recon, T_max, "Error-L2-" * problem.file_name * "-SLMsR-EdgeEvolved", "L2")
            Vis.writeError2file(errMSR_recon_H1, T_max, "Error-H1-" * problem.file_name * "-SLMsR-EdgeEvolved", "H1")

            Vis.writeSolution_all(solMSR_mapped, meshREF, problem.file_name * "-SLMsR-EdgeEvolved")
        else
            errMSR = PostProcess.error_L2(solREF, solMSR_mapped)[:]
            errMSR_H1 = PostProcess.error_H1(solREF, solMSR_mapped)[:]
            Vis.writeError2file(errMSR, T_max, "Error-L2-" * problem.file_name * "-SLMsR", "L2")
            Vis.writeError2file(errMSR_H1, T_max, "Error-H1-" * problem.file_name * "-SLMsR", "H1")

            Vis.writeSolution_all(solMSR_mapped, meshREF, problem.file_name * "-SLMsR-")
        end
        # ---------------------
        
        # ---------------------
        # Write a nodal basis function
        if true
            ind_basis = 35
            if reconstruct_edge
                Vis.writeNodalBasis_all_steps(solMSR, meshMSR, ind_basis, "Basis---node-" * string(ind_basis) * "-EdgeEvolved-" * problem.file_name)
            else
                Vis.writeBasis_all_steps(solMSR, meshMSR, ind_basis, "Basis---node-" * string(ind_basis) * "-" * problem.file_name)
            end
        end
        # ---------------------
end
# ---------------------------------------------------------------------------