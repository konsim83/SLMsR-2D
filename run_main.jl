if false
        include("reload_modules.jl")
end


# ---------------------------------------------------------------------------

import Parameter, Problem, FEM, Solver, PostProcess, Vis

# ---------------------------------------------------------------------------


compute_low = true
compute_ref = true
compute_ms = true

post_process = true


# ---------------------------------------------------------------------------
# -------   Problem Parameters   -------
T_max = 0.5


problem = Problem.Gaussian(T_max)

problem = Problem.Gaussian_1(T_max, 20)

problem = Problem.Gaussian_2(T_max, 1)
problem = Problem.Gaussian_2a(T_max, 1)

# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# -------   Mesh parameters   -------
n_edge_per_seg = 5
n_refinement = 3
n_edge_per_seg_f = 0


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 3

n_order_FEM_f = 1
n_order_quad_f = n_order_quad


time_step_method = 1

dt = 1/500
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
end


if compute_ms
        # -------   Build parameter structure   -------
        par_MsFEM = Parameter.Parameter_MsFEM(problem.T,
                                              dt,
                                              n_edge_per_seg,
                                              n_refinement,
                                              n_edge_per_seg_f,
                                              n_order_FEM_f,
                                              n_order_quad_f,
                                              time_step_method)


        # -------   Call the solver   -------
        @time solution_ms, mesh_collection = Solver.solve_MsFEM_periodic_square(par_MsFEM, problem)
end
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
if post_process
        solution_FEM_low_mapped = PostProcess.map_solution(solution_FEM_low, mesh_FEM_low, mesh_FEM_ref)
        solution_ms_mapped = PostProcess.map_solution(solution_ms, mesh_collection, mesh_FEM_ref)

        err_low = PostProcess.error_L2(solution_FEM_ref,
                                        solution_FEM_low_mapped)[:]
        err_ms = PostProcess.error_L2(solution_FEM_ref,
                                        solution_ms_mapped)[:]

        Vis.writeSolution_all(solution_FEM_ref, 
                                mesh_FEM_ref,
                                problem.file_name * "-FEM-ref")
        Vis.writeSolution_all(solution_FEM_low_mapped, 
                                mesh_FEM_ref, 
                                problem.file_name * "-FEM-low-mapped")
        Vis.writeSolution_all(solution_ms_mapped, 
                                mesh_FEM_ref, 
                                problem.file_name * "-MsFEM-mapped")
end
# ---------------------------------------------------------------------------