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


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 2

time_step_method = 1

dt = 1/500


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


if true
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
# ---------------------------------------------------------------------------