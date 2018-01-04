#include("Mesh.jl")
using PyPlot

P1 = [0.0 0.0 ; 1.0 0.0 ; 1.0 1.0];
ms1 = Mesh.mesh_triangle_uniform_edges(P1, 3);

P2 = [0.0 0.0 ;1.0 1.0 ; 0 1.0]+1;
ms2 = Mesh.mesh_triangle_uniform_edges(P2, 3);



function add_meshes(mesh1 :: Mesh.TriMesh, mesh2 :: Mesh.TriMesh)

    p1, c1 = copy(mesh1.point), copy(mesh1.cell)
    p2, c2 = copy(mesh2.point), copy(mesh2.cell)

    n_p1 = mesh1.n_point
    n_p2 = mesh2.n_point
   
    c2 += n_p1

    P = [p1 ; p2]
    C = [c1 ; c2]
    
    #n_found = 0
    ind_delete = []
    for j in 1:n_p1
        #p1_in_p2 = (Bool[ p1[j,:] == p2[i,:] for i in 1:n_p2 ])
        p1_in_p2 = (Bool[ sum(abs.(p1[j,:] - p2[i,:]))<1e-10 for i in 1:n_p2 ])
        display(p1_in_p2)
        # This index in c2 must be set to j
        ind_p1_in_p2 = find(p1_in_p2)
        
        length(ind_p1_in_p2)>1 ? error("something is wrong...") :
            if ~isempty(ind_p1_in_p2)
                push!(ind_delete, ind_p1_in_p2[1])
                
                c2[c2.==ind_p1_in_p2+n_p1] = -j # set it to -j to easily identify changed points
                # increase by one if a point is found
                #n_found += sum(p1_in_p2)
            end
    end
    
    p2 = p2[setdiff(1:n_p2,ind_delete),:]
    #c2[c2.>0] = c2[c2.>0]  + n_p1 #- n_found
    c2 = abs.(c2)

    P = [p1 ; p2]
    C = [c1 ; c2]

    C = unique(C,1)
       
    return P, C
end


function plot_mesh(P :: Array{Float64,2}, C :: Array{Int64,2})

    fig = matplotlib[:pyplot][:figure]("2D Mesh Plot", figsize = (10,10))
    
    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")
    
    # Connectivity list -1 for Python
    tri = ax[:triplot](P[:,1], P[:,2], C - 1 )
    setp(tri, linestyle = "-",
         marker = "None",
         linewidth = 1)
    
    fig[:canvas][:draw]()
    
    return nothing
end





#=
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
@time A = FEM.assemble_advection(m, d, r, q, par, p, 1)
#@time D = FEM.assemble_diffusion(m, d, r, q, par, p, 1)
println("---------------------------\n\n")



# -----------------------------------------------------------------------------------------------------------------
ind_cell = FEM.get_dof_elem(d, m, 1:d.n_elem)
ind = vec(ind_cell[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
ind_test = vec(transpose(repmat(ind_cell, 1, size(ind_cell,2))))
ind_lin = sub2ind((d.n_true_dof,d.n_true_dof), ind_test, ind)

#M1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
A1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#D1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
# -----------------------------------------------------------------------------------------------------------------

println("In place assembly...")
#@time FEM.assemble_mass!(M1, m, d, r, q, par, p)
@time FEM.assemble_advection!(A1, m, d, r, q, par, p, 1)
#@time FEM.assemble_diffusion!(D1, m, d, r, q, par, p, 1)
println("---------------------------\n\n")
=#
