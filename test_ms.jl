if true
    include("Parameter.jl")
    include("Problem.jl")
    include("Quad.jl")
    include("Mesh.jl")
    include("FEM.jl")
    include("Time_integrator.jl")
    include("Solver.jl")
end

n_edge_per_seg = 10
n_edge_per_seg_f = 10
n_order_quad_f = 4
n_order_FEM_f = 1


par_ms = Parameter.Parameter_MsFEM(1.0,
                                   1/500,
                                   n_edge_per_seg,
                                   n_edge_per_seg_f,
                                   n_order_FEM_f,
                                   n_order_quad_f,
                                   1)


m_coarse = Mesh.mesh_unit_square(n_edge_per_seg)
m_collection = Mesh.TriMesh_collection(m_coarse, n_edge_per_seg_f)
#=
mesh_simplex = Mesh.mesh_unit_simplex_uniform_edges(par.n_edge_per_seg)
mesh_collection = Mesh.Trimesh_collection(m_coarse, mesh_simplex)
=#


r = FEM.RefEl_Pk{1}()
r_f = FEM.RefEl_Pk{n_order_FEM_f}()


q_f = Quad.Quad_simplex(n_order_quad_f)


d_collection = FEM.Dof_collection{r_f.n_order}(m_collection)


s = FEM.Solution_MsFEM(d_collection, par_ms)


p = Problem.Gaussian(1.0);
p_f = Problem.BasisFun(1, 1.0);


i = 1;
tstep = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(d_collection.dof_f[1],
                                                                             m_collection.mesh_f[i],
                                                                             p_f)
# ===============================================
