#=
include("Parameter.jl")
include("Problem.jl")
include("Quad.jl")
include("Mesh.jl")
include("FEM.jl")
include("Time_integrator.jl")
include("Solver.jl")
=#

n_edge_per_seg = 10
n_order_quad = 4
n_order_FEM = 1


println("Start timing problem setup...")
begin
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

# ===============================================




println("Standard assembly...")
#@time M = FEM.assemble_mass(m, d, r, q, par, p, 1)
@time A = FEM.assemble_advection(m, d, r, q, par, p, 1)
#@time D = FEM.assemble_diffusion(m, d, r, q, par, p, 1)
println("---------------------------\n\n")


i = FEM.get_dof_elem(d, m, 1:d.n_elem)
ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
ind_test = vec(transpose(repmat(i, 1, size(i,2))))
lin = sub2ind((d.n_true_dof,d.n_true_dof), ind_test, ind)


#M1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
A1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#D1 = sparse(ind_test, ind, zeros(Float64, length(ind)), d.n_true_dof, d.n_true_dof)
#println("\n")

println("In place assembly...")
#@time FEM.assemble_mass!(M1, m, d, r, q, par, p, 1)
@time FEM.assemble_advection!(A1, m, d, r, q, par, p, 1)
#@time FEM.assemble_diffusion!(D1, m, d, r, q, par, p, 1)
println("---------------------------\n\n")



# ===============================================
function modify_ansatzfunction_v(v :: Array{Float64,2}, DF :: Array{Float64,3}, DPhi :: Array{Float64,3})

    DPhi_mod = Array{Float64,2}(size(DPhi,1), size(DPhi,2))

    for i=1:size(DPhi,2)
        DPhi_mod[:,i] = DPhi[:,i,:] * DF[i,:,:] * v[i,:]
    end
    
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


# ===============================================
function mesh_expand(DPhi :: Array{Float64,3}, DF :: Array{Float64,4})

    DPhi_mesh = Array{Float64,4}(size(DPhi,1),size(DPhi,2), size(DPhi,3), size(DF,4))

    for i in 1:size(DPhi,2)
        for k in 1:size(DF,4)
            DPhi_mesh[:,i,:,k] = DPhi[:,i,:] * DF[i,:,:,k]
        end
    end

    return DPhi_mesh
end

function modify_ansatzfunction_v(v :: Array{Float64,2}, DPhi_mesh :: Array{Float64,3})

    DPhi_mod = Array{Float64,2}(size(DPhi_mesh,1),size(DPhi_mesh,2))

    for i=1:size(DPhi_mesh,1)
        DPhi_mod[i,:] = sum(DPhi_mesh[i,:,:] .* v,2)
    end
    
    return DPhi_mod
end

# ===============================================


# ----------------
println("Starting assembly test....")
i = FEM.get_dof_elem(d, m, 1:d.n_elem)
ind_test = vec(transpose(repmat(i, 1, size(i,2))))
ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
lin = sub2ind((d.n_true_dof,d.n_true_dof), ind_test, ind)

n = r.n_node


# fixed quantities for mesh
weight_elem = Mesh.map_ref_point_grad_det(m, q.point, 1:d.n_elem)
DF = Mesh.map_ref_point_grad_inv(m, q.point, 1:d.n_elem);
x = Mesh.map_ref_point(m, q.point, 1:d.n_elem)



# ---   Advection   ---------------------------------------
#-------------------------------------------------------------
println("... Advection routines:")

velocity = Problem.velocity(p, 1.0, x)

DPhi = FEM.eval_grad(r, q.point)
Phi_test = FEM.eval(r, q.point)*diagm(q.weight)

DPhi_mesh = mesh_expand(DPhi, DF)


DPhi_mod = 0
DPhi_mod_new = 0

@time begin
    a1_loc = Array{Float64,3}(n,n,m.n_cell)
    
    #println("****** First loop *******")
    for k = 1:m.n_cell
        #DPhi_mod = modify_ansatzfunction_v(velocity[:,:,k], DF[:,:,:,k], DPhi)
        DPhi_mod = modify_ansatzfunction_v(velocity[:,:,k], DPhi_mesh[:,:,:,k])
        a1_loc[:,:,k] = Phi_test * diagm(weight_elem[:,k])  * DPhi_mod'
    end

    A1 = sparse(ind_test, ind, vec(a1_loc), d.n_true_dof, d.n_true_dof)
end


A2 = sparse(d.ind_test, d.ind, zeros(Float64, length(d.ind)), d.n_true_dof, d.n_true_dof)
@time begin
    a2_loc = Array{Float64,2}(n,n)
    #println("****** Second loop *******")
    for k = 1:m.n_cell
        #DPhi_mod = modify_ansatzfunction_v(velocity[:,:,k], DF[:,:,:,k], DPhi)
        DPhi_mod = modify_ansatzfunction_v(velocity[:,:,k], DPhi_mesh[:,:,:,k])
        a2_loc[:,:] = Phi_test * diagm(weight_elem[:,k]) * DPhi_mod'
        for l in 1:9
            A2[d.ind_test[9*(k-1)+l],d.ind[9*(k-1)+l]] += a2_loc[l]
        end
    end
end


A3 = sparse(d.ind_test, d.ind, zeros(Float64, length(d.ind)), d.n_true_dof, d.n_true_dof)
@time begin
    a3_loc = Array{Float64,2}(n,n)
    #println("****** Second loop *******")
    for k = 1:m.n_cell
        for j = 1:length(q.weight)

            
            for ii_test in 1:n
                for ii in 1:n
                    a3_loc[ii_test,ii] = Phi_test[ii_test,j] * weight_elem[j,k] * (DPhi[ii,j,1]*velocity[j,1,k] + DPhi[ii,j,2]*velocity[j,2,k])
                end
            end

            
            #a3_loc[:,:] = Phi_test[:,j] * weight_elem[j,k] * (DPhi[:,j,:]*velocity[j,:,k])'
            for l in 1:9
                A3[d.ind_test[9*(k-1)+l],d.ind[9*(k-1)+l]] += a3_loc[l]
            end
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
