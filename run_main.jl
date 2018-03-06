if true
        include("Parameter.jl")
        include("Geometry.jl")
        include("Problem.jl")
        include("Quad.jl")
        include("Mesh.jl")
        include("FiniteDiff.jl")
        include("FEM.jl")
        include("Time_integrator.jl")
        include("Solver.jl")
        include("PostProcess.jl")
        include("Vis.jl")
end
# ---------------------------------------------------------------------------


compute_low = true
compute_high = true
compute_ms = true

post_process = true


# ---------------------------------------------------------------------------
# -------   Problem Parameters   -------
T_max = 0.5

problem = Problem.Gaussian(T_max)

# problem = Problem.Gaussian_1(T_max)

# problem = Problem.Gaussian(T_max, 1)
# problem = Problem.Gaussiana(T_max, 1)

# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# -------   Mesh parameters   -------
n_edge_per_seg = 5
n_refinement = 3
n_edge_per_seg_f = 0


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 2

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

if compute_high
        # -------   Build parameter structure   -------
        par_FEM_high = Parameter.Parameter_FEM(problem.T,
                                                 dt,
                                                 n_edge_per_seg,
                                                 n_refinement,
                                                 n_order_FEM,
                                                 n_order_quad,
                                                 time_step_method)

        # -------   Call the solver   -------
        @time solution_FEM_high, mesh_FEM_high = Solver.solve_FEM_periodic_square(par_FEM_high, problem)
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
        solution_FEM_low_mapped = PostProcess.map_solution(solution_FEM_low, mesh_FEM_low, mesh_FEM_high)
        solution_ms_mapped = PostProcess.map_solution(solution_ms, mesh_collection, mesh_FEM_high)

        # err_low_2 = PostProcess.error_L2(solution_FEM_high, solution_FEM_low_mapped)[:]
        # err_ms_2 = PostProcess.error_L2(solution_FEM_high, solution_ms_mapped)[:]
        # Vis.writeSolution_all(solution_FEM_high, mesh_FEM_high, "Gaussian-FEM-high---classic")
        # Vis.writeSolution_all(solution_FEM_low_mapped, mesh_FEM_high, "Gaussian-FEM-low-mapped---classic")
        # Vis.writeSolution_all(solution_ms_mapped, mesh_FEM_high, "Gaussian-MsFEM-mapped---classic")

        err_low_2a = PostProcess.error_L2(solution_FEM_high, solution_FEM_low_mapped)[:]
        err_ms_2a = PostProcess.error_L2(solution_FEM_high, solution_ms_mapped)[:]
        Vis.writeSolution_all(solution_FEM_high, mesh_FEM_high, "Gaussian-FEM-high---stream")
        Vis.writeSolution_all(solution_FEM_low_mapped, mesh_FEM_high, "Gaussian-FEM-low-mapped---stream")
        Vis.writeSolution_all(solution_ms_mapped, mesh_FEM_high, "Gaussian-MsFEM-mapped---stream")

end
---------------------------------------------------------------------------