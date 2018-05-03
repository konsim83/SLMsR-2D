# ---------------------------------------------------------------------------

import Parameter, Problem, FEM, Solver, PostProcess, Vis, Reconstruction

# ---------------------------------------------------------------------------


compute_low = true
compute_ref = true
compute_ms = false
compute_ms_reconstruction = false

post_process = false


# ---------------------------------------------------------------------------
# -------   Problem Parameters   -------
T_max = 0.5


# problem = Problem.Gaussian(T_max)

# problem = Problem.Gaussian_1(T_max, 20)

# problem = Problem.Gaussian_2(T_max, 1)
# problem = Problem.Gaussian_2a(T_max, 1)



# Multiscale diffusion, constant advection
# problem = Problem.Gaussian_R_1(T_max, 20)


# Constant low diffusion, solenoidal, traveling vortex
# problem = Problem.Gaussian_R_2(T_max, 0.3 , 20)
# problem = Problem.Gaussian_R_2_conserv(T_max, 0.3 , 30)


# Multiscale diffusion, solenoidal, traveling vortex
# problem = Problem.Gaussian_R_3(T_max, 0.1 , 30)
# problem = Problem.Gaussian_R_3_conserv(T_max, 0.1 , 30)


# Multiscale diffusion, divergent, traveling vortex
problem = Problem.Gaussian_R_4(T_max, 0.05 , 30)
# problem = Problem.Gaussian_R_4_conserv(T_max, 0.05 , 30)

# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# -------   Mesh parameters   -------
n_edge_per_seg = 4
n_refinement = 4
n_edge_per_seg_f = 0


# -------   FEM parameters   -------
n_order_FEM = 1
n_order_quad = 3

n_order_FEM_f = 1
n_order_quad_f = n_order_quad


time_step_method = 1

dt = 1/250
n_steps_f = 5

k = [1.0 ; 1.0 ; 1.0] * 0.0001
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
                                                n_steps_f,
                                                n_edge_per_seg,
                                                n_refinement,
                                                n_edge_per_seg_f,
                                                n_order_FEM_f,
                                                n_order_quad_f,
                                                time_step_method,
                                                k)


        # -------   Call the solver   -------
        @time solution_ms, mesh_collection = Solver.solve_MsFEM_periodic_square(par_MsFEM, problem)
end


if compute_ms_reconstruction
        # -------   Build parameter structure   -------
        par_MsFEM_r = Parameter.Parameter_MsFEM(problem.T,
                                                dt,
                                                n_steps_f,
                                                n_edge_per_seg,
                                                n_refinement,
                                                n_edge_per_seg_f,
                                                n_order_FEM_f,
                                                n_order_quad_f,
                                                time_step_method,
                                                k)


        # -------   Call the solver   -------
        @time solution_ms_r, mesh_collection_r = Solver.solve_MsFEM_periodic_square_reconstruction(par_MsFEM_r, problem)
end
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
if post_process
        solution_FEM_low_mapped = PostProcess.map_solution(solution_FEM_low, mesh_FEM_low, mesh_FEM_ref)
        # solution_ms_mapped = PostProcess.map_solution(solution_ms, mesh_collection, mesh_FEM_ref)
        solution_ms_r_mapped = PostProcess.map_solution(solution_ms_r, mesh_collection_r, mesh_FEM_ref)

        err_low = PostProcess.error_L2(solution_FEM_ref,
                                        solution_FEM_low_mapped)[:]
        # err_ms = PostProcess.error_L2(solution_FEM_ref,
        #                                 solution_ms_mapped)[:]
        err_ms_r = PostProcess.error_L2(solution_FEM_ref,
                                        solution_ms_r_mapped)[:]

        Vis.writeSolution_all(solution_FEM_ref, 
                                mesh_FEM_ref,
                                problem.file_name * "-FEM-ref")
        Vis.writeSolution_all(solution_FEM_low_mapped, 
                                mesh_FEM_ref, 
                                problem.file_name * "-FEM-low-mapped")
        # Vis.writeSolution_all(solution_ms_mapped, 
        #                         mesh_FEM_ref, 
        #                         problem.file_name * "-MsFEM-mapped")
        Vis.writeSolution_all(solution_ms_r_mapped, 
                                mesh_FEM_ref, 
                                problem.file_name * "-MsFEM_r-mapped")
end
# ---------------------------------------------------------------------------



#=
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
par = Parameter.Parameter_MsFEM(problem.T,
                                dt,
                                n_edge_per_seg,
                                n_refinement,
                                n_edge_per_seg_f,
                                n_order_FEM_f,
                                n_order_quad_f,
                                time_step_method,
                                k)

# Build mesh of unit square (0,1)x(0,1)
m_coarse = Mesh.mesh_unit_square(par.n_edge_per_seg)
m_simplex = Mesh.mesh_unit_simplex()
if par.n_refinement>0
# If this is not the case then the conformal MsFEM is just the
# standard FEM
m_simplex = Mesh.refine_rg(m_simplex, par.n_refinement)

# Subdividing edges messes up boundary markers. We need to correct
# that.
ind_point_boundary = sort(unique(m_simplex.edge[:,m_simplex.edge_marker.!=0]))
m_simplex.point_marker[:] = zeros(Int, size(m_simplex.point_marker))
m_simplex.point_marker[ind_point_boundary] = ones(Int, size(ind_point_boundary))
end

mesh_collection = Mesh.TriMesh_collection(m_coarse, m_simplex)
#mesh_collection = Mesh.TriMesh_collection(m_coarse, par.n_edge_per_seg_f)


# Set up reference element
ref_el = FEM.RefEl_Pk(1)
ref_el_f = FEM.RefEl_Pk(par.n_order_FEM_f)

# Set up quadrature rule
quad_f = Quad.Quad_simplex(par.n_order_quad_f)

# Set up degrees of freedom handler
periodicityInfo = Mesh.DoublePeriodicUnitSquare()
dof = FEM.Dof_Pk_periodic(mesh_collection.mesh, problem, periodicityInfo, 1)
problem_f = Array{Problem.AbstractBasisProblem,1}(mesh_collection.mesh.n_cell)
for i_cell in 1:mesh_collection.mesh.n_cell
# Set up local problem by geometry
tri = Geometry.Triangle(FEM.map_ref_point(dof, ref_el.node, i_cell))
problem_f[i_cell] = Problem.BasisFun(problem, tri)
end
dof_collection = FEM.Dof_collection(mesh_collection, dof, problem_f, ref_el_f.n_order)

# Set up solution structure
solution = FEM.Solution_MsFEM(dof_collection, par)


# ----------------------------------------------------------------
# ----------------------------------------------------------------
# Setup initial data for basis functions
timeStepper_local = Array{TimeIntegrator.ImplEuler,1}(mesh_collection.mesh.n_cell)
for ind_cell in 1:mesh_collection.mesh.n_cell
u_init_tmp = Problem.u_init(problem_f[ind_cell], mesh_collection.mesh_f[ind_cell].point)
solution.phi_1[ind_cell][:,1] = u_init_tmp[:,1]
solution.phi_2[ind_cell][:,1] = u_init_tmp[:,2]
solution.phi_3[ind_cell][:,1] = u_init_tmp[:,3]

timeStepper_local[ind_cell] = TimeIntegrator.ImplEuler(dof_collection.dof_f[ind_cell],
                                                        mesh_collection.mesh_f[ind_cell],
                                                        problem_f[ind_cell])
end

# Setup initial data
solution.u[1:mesh_collection.mesh.n_point,1] = Problem.u_init(problem, mesh_collection.mesh.point)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

timeStepper = TimeIntegrator.ImplEuler(dof_collection.dof,
                                            mesh_collection.mesh,
                                            problem)

for k_time=1:1#par.n_steps
        # Time at index k
        Time = (k_time-1)*par.dt

        Reconstruction.SemiLagrange_L2_opt!(solution, 
                                                timeStepper_local,
                                                mesh_collection,
                                                dof_collection,
                                                ref_el_f,
                                                quad_f,
                                                par,
                                                problem,
                                                problem_f,
                                                Time + par.dt, k_time + 1)

        M, Mt = FEM.assemble_mass(solution,
                                   mesh_collection,
                                   dof_collection,
                                   ref_el, ref_el_f,
                                   quad_f,
                                   problem,
                                   k_time+1)

        A = FEM.assemble_advection(solution,
                                   mesh_collection,
                                   dof_collection,
                                   ref_el, ref_el_f,
                                   quad_f,
                                   problem,
                                   k_time+1, Time + par.dt)

        D = FEM.assemble_diffusion(solution,
                                   mesh_collection,
                                   dof_collection,
                                   ref_el, ref_el_f,
                                   quad_f,
                                   problem,
                                   k_time+1, Time + par.dt)

        # Zero forcing
        f = zeros(mesh_collection.mesh.n_point)

        # Set the system matrices 
        TimeIntegrator.updateSystem!(timeStepper.systemData, M, D-A-Mt, f, dof, solution.u[:,k_time], par.dt)

        # Make a single time step
        TimeIntegrator.makeStep!(timeStepper, dof, view(solution.u, :, k_time+1), solution.u[:,k_time])
end
=#