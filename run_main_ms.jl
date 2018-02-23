if false
	include("Parameter.jl")
	include("Geometry.jl")
	include("Problem.jl")
	include("Quad.jl")
	include("Mesh.jl")
	include("FEM.jl")
	include("Time_integrator.jl")
	include("Solver.jl")
	include("FiniteDiff.jl")
	include("PostProcess.jl")
	include("Vis.jl")
end
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
"""
Initialize the problem.
"""

# -------   Problem Parameters   -------
T_max = 0.5

problem = Problem.Gaussian(T_max)

# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
"""
Coarse standard FEM. Solves Advection diffusion equation on coarse mesh.
"""

# -------   Mesh parameters   -------
n_edge_per_seg = 5
n_refinement = 4
n_edge_per_seg_f = 0


# -------   FEM parameters   -------
n_order_FEM_f = 1
n_order_quad_f = 2

time_step_method = 1

dt = 1/500


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
# ---------------------------------------------------------------------------
