include("Parameter.jl")
include("Problem.jl")
include("Quad.jl")
include("Mesh.jl")
include("FEM.jl")
include("Time_integrator.jl")
include("Solver.jl")


n_edge_per_seg = 2^2
n_order_quad = 1
n_order_FEM = 1

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
    q = Quad.Quad_simplex(3);
    
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
# ===============================================




println("Standard assembly...")
#@time M = FEM.assemble_mass(m, d, r, q, par, p, 1)
#@time A = FEM.assemble_advection(m, d, r, q, par, p, 1)
#@time D = FEM.assemble_diffusion(m, d, r, q, par, p, 1)
println("---------------------------\n\n")


i = FEM.get_dof_elem(d, m, 1:d.n_elem)
ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
ind_test = vec(transpose(repmat(i, 1, size(i,2))))
lin = sub2ind((d.n_true_dof,d.n_true_dof), ind_test, ind)


M1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#A1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#D1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#println("\n")

println("In place assembly...")
#@time FEM.assemble_mass!(M1, m, d, r, q, par, p, 1)
#@time FEM.assemble_advection!(A1, sol, m, d, r, q, par, p, 1)
#@time FEM.assemble_diffusion!(D1, sol, m, d, r, q, par, p, 1)
println("---------------------------\n\n")



# ===============================================
function modify_ansatzfunction_v(v :: Array{Float64,2}, DF :: Array{Float64,3}, DPhi :: Array{Float64,3})

    DPhi_mod = Array{Float64,2}(size(DPhi,1), size(DPhi,2))

    for i=1:size(DPhi,1)
        DPhi_mod[i,:] = DPhi[i,:,:] * DF[i,:,:] * v[i,:]
    end
    
    return DPhi_mod
end


function modify_ansatzfunction_v_new(v :: Array{Float64,2}, DF :: Array{Float64,3}, DPhi :: Array{Float64,3})

    DPhi_mod = Array{Float64,2}(size(DPhi,1), size(DPhi,2))

    D = reshape([DPhi[:,1,:] * DF[1,:,:] ; DPhi[:,2,:] * DF[2,:,:] ; DPhi[:,3,:] * DF[3,:,:]], size(DPhi,1), size(DPhi,2), 2)
    DPhi_mod = broadcast(*, v[:,1]', D[:,:,1]) + broadcast(*, v[:,2]', D[:,:,2])

    return DPhi_mod
end
# ===============================================


# ===============================================
function modify_ansatzfunction_d(d :: Array{Float64,3}, DF :: Array{Float64,3},
                                 Phi_test :: Array{Float64,3}, Phi :: Array{Float64,3},
                                 q_weight :: Array{Float64,1}, weight_elem :: Array{Float64,1})
    
    Phi_mod = zeros(size(Phi,1), size(Phi,2))
    
    for i=1:size(Phi,1)
        Phi_mod += (Phi_test[i,:,:] * DF[i,:,:]) * q_weight[i]*weight_elem[i] *  (Phi[i,:,:] * DF[i,:,:] * d[i,:,:])'
    end
    
    return -Phi_mod
end
# ===============================================


# ----------------
println("Starting assembly test....")
i = FEM.get_dof_elem(d, m, 1:d.n_elem)
ind_test = vec(transpose(repmat(i, 1, size(i,2))))
ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')

n = r.n_node
mat_local = zeros(n, n, m.n_cell)

@time x = Mesh.map_ref_point(m, q.point, 1:d.n_elem)



# ---   Advection   ---------------------------------------
#-------------------------------------------------------------
println("... Advection routines:")
@time velocity = Problem.velocity(p, 1.0, x)

@time DF = Mesh.map_ref_point_grad_inv(m, q.point, 1:d.n_elem);
@time weight_elem = Mesh.map_ref_point_grad_det(m, q.point, 1:d.n_elem)

DPhi = FEM.eval_grad(r, q.point)
Phi_test = diagm(q.weight) * FEM.eval(r, q.point)

DPhi_mod = 0
DPhi_mod_new = 0

begin
    for k = 1#:m.n_cell
        #k=1
        @time DPhi_mod = modify_ansatzfunction_v(velocity[:,:,k], DF[:,:,:,k], DPhi)
        @time mat_local[:,:,k] = Phi_test' * diagm(weight_elem[:,k]) * DPhi_mod
    end

    A1 = sparse(ind_test, ind, vec(mat_local), d.n_true_dof, d.n_true_dof)
end

A2 = sparse(d.ind_test, d.ind, zeros(Float64, length(d.ind)), d.n_true_dof, d.n_true_dof)
begin
    mat_loc = Array{Float64,2}(3,3)
    
    for k = 1#:m.n_cell
        @time DPhi_mod_new = modify_ansatzfunction_v_new(velocity[:,:,k], DF[:,:,:,k], DPhi)
        @time mat_loc[:,:] = Phi_test' * diagm(weight_elem[:,k]) * DPhi_mod_new
        for l in 1:9
            A2[d.ind_test[9*(k-1)+l],d.ind[9*(k-1)+l]] += mat_loc[l]
        end
    end
end

println("\n")
#-------------------------------------------------------------





#=
# ---   Diffusion   ----------------------------------------
#-------------------------------------------------------------
println("... Diffusion routines:")
@time diff = Problem.diffusion(p, 1.0, x)


DPhi = FEM.eval(r, q.point);
DPhi_test = FEM.eval(r, q.point);

@time DF = Mesh.map_ref_point_grad_inv(m, q.point, 1:d.n_elem);
@time weight_elem = Mesh.map_ref_point_grad_det(m, q.point, 1:d.n_elem)

@time begin
    for k = 1:m.n_cell    
        mat_local[:,:,k] = modify_ansatzfunction_d(diff[:,:,:,k], DF[:,:,:,k], DPhi_test, DPhi, q.weight, weight_elem[:,k])
    end
end
println("\n")
#-------------------------------------------------------------
=#

#D1 = sparse(ind_test, ind, vec(mat_local), d.n_true_dof, d.n_true_dof)

println("---------------------------\n\n")
