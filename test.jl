include("reload.jl");

m = Mesh.mesh_unit_square(2);
r = FEM.RefEl_Pk{1}();
d = FEM.Dof_Pk_periodic_square{r.n_order}(m);
p = Problem.Gaussian(1.0);
q = Quad.Quad_simplex(3);


# ===============================================

function modify_ansatzfunction_v(v :: Array{Float64,2}, DF :: Array{Float64,3}, Phi :: Array{Float64,3})

    Phi_mod = Array{Float64,2}(size(Phi,1), size(Phi,2))

    for i=1:size(Phi,1)
        Phi_mod[i,:] = Phi[i,:,:] * DF[i,:,:] * v[i,:]
    end
    
    return Phi_mod
end



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
i = FEM.get_dof_elem(d, m, 1:d.n_elem)
ind_test = vec(transpose(repmat(i, 1, size(i,2))))
ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')



n = r.n_node
mat_local = zeros(n, n, m.n_cell)
    
x = Mesh.map_ref_point(m, q.point, 1:d.n_elem);
diff = Problem.diffusion(p, 1.0, x)


Phi = FEM.eval(r, q.point);
Phi_test = FEM.eval(r, q.point);

DF = Mesh.map_ref_point_grad_inv(m, q.point, 1:d.n_elem);
weight_elem = Mesh.map_ref_point_grad_det(m, q.point, 1:d.n_elem)

for k = 1:m.n_cell    
    mat_local[:,:,k] = modify_ansatzfunction_d(diff[:,:,:,k], DF[:,:,:,k], Phi_test, Phi, q.weight, weight_elem[:,k])
end


Mat_global = sparse(ind_test, ind, vec(mat_local), d.n_true_dof, d.n_true_dof)

