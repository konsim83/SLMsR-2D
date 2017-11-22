"""
Coarse standard FEM. Solves Advection diffusion equation on coarse mesh.
"""

# -------   Mesh parameters   -------
n_edge_per_seg = 5


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 2

time_step_method = 1

dt = 1/500



# -------   Build parameter structure   -------
par_FEM_coarse = Parameter.Parameter_FEM(problem.T,
                                         dt,
                                         n_edge_per_seg,
                                         n_n_order_FEM,
                                         n_order_quad,
                                         time_step_method)


# -------   Call the solver   -------
@time solution, mesh = FEM.solve_FEM(par_FEM_coarse, problem)
