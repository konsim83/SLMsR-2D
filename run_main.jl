include("Parameter.jl")
include("Geometry.jl")
include("Problem.jl")
include("Quad.jl")
include("Mesh.jl")
include("FiniteDiff.jl")
include("FEM.jl")
include("Time_integrator.jl")
include("Solver.jl")
include("Vis.jl")
# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
"""
Initialize the problem.
"""

# -------   Problem Parameters   -------
T_max = 1.0

problem = Problem.Gaussian(T_max)

# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
"""
Coarse standard FEM. Solves Advection diffusion equation on coarse mesh.
"""

# -------   Mesh parameters   -------
n_edge_per_seg = 10
n_refinement = 1


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 2

time_step_method = 1

dt = 1/500


# -------   Build parameter structure   -------
par_FEM_coarse = Parameter.Parameter_FEM(problem.T,
                                         dt,
                                         n_edge_per_seg,
                                         n_refinement,
                                         n_order_FEM,
                                         n_order_quad,
                                         time_step_method)


# -------   Call the solver   -------
@time solution_FEM, mesh_FEM = Solver.solve_FEM_periodic_square(par_FEM_coarse, problem)
# ---------------------------------------------------------------------------
