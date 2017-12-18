include("Parameter.jl")
include("Problem.jl")
include("Quad.jl")
include("Mesh.jl")
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

ind_basis = 1
problem = Problem.BasisFun(ind_basis, T_max)

# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
"""
Standard FEM. Solves Advection diffusion equation for the basis functions.
"""

# -------   Mesh parameters   -------
n_edge_per_seg = 20
#n_edge_per_seg = 50


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 2

time_step_method = 1

dt = 1/500


# -------   Build parameter structure   -------
par_FEM_coarse = Parameter.Parameter_FEM(problem.T,
                                         dt,
                                         n_edge_per_seg,
                                         n_order_FEM,
                                         n_order_quad,
                                         time_step_method)


# -------   Call the solver   -------
@time solution, mesh = Solver.solve_FEM_simplex(par_FEM_coarse, problem)
# ---------------------------------------------------------------------------