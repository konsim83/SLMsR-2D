if true
    include("Parameter.jl")
    include("Problem.jl")
    include("Quad.jl")
    include("Mesh.jl")
    include("FEM.jl")
    include("Time_integrator.jl")
    include("Solver.jl")
end

n_edge_per_seg = 5
n_order_quad = 4
n_order_FEM = 1


if true
    println("Start timing problem setup...")
    @time begin
        println("\tMeshing...")
        m = Mesh.mesh_unit_square(n_edge_per_seg);
        
        println("\tReference element...")
        r = FEM.RefEl_Pk{n_order_FEM}();
        
        println("\tDof handler...")
        d = FEM.Dof_Pk_periodic_square{n_order_FEM}(m);
        
        println("\tProblem structure...")
        p = Problem.Gaussian(1.0);
        
        println("\tQuadrature...")
        q = Quad.Quad_simplex(n_order_quad);
        
        par = Parameter.Parameter_FEM(1.0,
                                      1.0/1000,
                                      n_edge_per_seg,
                                      n_order_FEM,
                                      n_order_quad,
                                      1)
        
        #sol = FEM.Solution_FEM(d, par)
        
        println("\tTime integrator...")
        tstep = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(d, m, p)
    end
    println("----------------------------\n")
end
# ===============================================




println("Standard assembly...")
#@time M = FEM.assemble_mass(m, d, r, q, par, p)
#@time A = FEM.assemble_advection(m, d, r, q, par, p, 1)
@time D = FEM.assemble_diffusion(m, d, r, q, par, p, 1)
println("---------------------------\n\n")




#M1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#A1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
D1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#println("\n")

println("In place assembly...")
#@time FEM.assemble_mass!(M1, m, d, r, q, par, p)
#@time FEM.assemble_advection!(A1, m, d, r, q, par, p, 1)
@time FEM.assemble_diffusion!(D1, m, d, r, q, par, p, 1)
println("---------------------------\n\n")
