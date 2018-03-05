include("Parameter.jl")
include("Geometry.jl")
include("Problem.jl")
include("Quad.jl")
include("Mesh.jl")
include("FiniteDiff.jl")
include("FEM.jl")


par_global = Parameter.Parameter_FEM(1.0,
		                                1/500,
		                                2,
		                                0,
		                                1,
		                                3,
		                                1)

m = Mesh.mesh_unit_square(10)
r = FEM.RefEl_Pk(1)
p = Problem.Gaussian(1.0)
per = Mesh.DoublePeriodicUnitSquare()
q = Quad.Quad_simplex(3)
d = FEM.Dof_Pk_periodic(m, p, per, 1)

# ------------------------
# -------   Mass   -------
# ------------------------
function assemble_mass(mesh :: Mesh.TriangleMesh.TriMesh,
                       dof :: FEM.AbstractDof,
                       ref_el :: FEM.RefEl_Pk{N},
                       quad :: Quad.AbstractQuad,
                       problem :: Problem.AbstractProblem) where N
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_m(mesh, ref_el, dof, quad)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)

    return Mat_global
end


function assemble_elem_m(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk{N},
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad) where N
    
    n = ref_el.n_node
    mass = Array{Float64,3}(n, n, mesh.n_cell)
    
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    Phi = FEM.shapeFun(ref_el, quad.point)
    Phi_test = FEM.shapeFun(ref_el, quad.point) * diagm(quad.weight)

    for k = 1:mesh.n_cell
        mass[:,:,k] = Phi_test * weight_elem[k] * Phi'
    end

    return mass
end # end function



# -----------------------------
# -------   Advection   -------
# -----------------------------
function assemble_advection(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.AbstractProblem,
                            k_time :: Int64)

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_a(mesh, ref_el, dof, quad, problem, time)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_a(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time_idx :: Float64)
    
    n = ref_el.n_node
    a = zeros(n, n, mesh.n_cell)
    
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    velocity = Problem.velocity(problem, time_idx, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem)
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    DPhi = FEM.shapeFun_grad(ref_el, quad.point)
    Phi_test = FEM.shapeFun(ref_el, quad.point)

    for k = 1:mesh.n_cell
    	for j=1:length(quad.weight)
    		for i=1:ref_el.n_node
        		a[i,j,k] = weight_elem[k] * quad.weight[j] * Phi_test[i,j]  * dot(velocity[k][j],(DF[k]'*DPhi[i,j]))
        	end
        end
    end

    return a
end # end function



# -----------------------------
# -------   Diffusion   -------
# -----------------------------
function assemble_diffusion(mesh :: Mesh.TriangleMesh.TriMesh,
                            dof :: FEM.AbstractDof,
                            ref_el :: FEM.RefEl_Pk,
                            quad :: Quad.AbstractQuad,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.AbstractProblem,
                            k_time :: Int64)

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_d(mesh, ref_el, dof, quad, problem, time)

    Mat_global = sparse(dof.ind_test, dof.ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    return Mat_global
end


function assemble_elem_d(mesh :: Mesh.TriangleMesh.TriMesh,
                         ref_el :: FEM.RefEl_Pk,
                         dof :: FEM.AbstractDof,
                         quad :: Quad.AbstractQuad,
                         problem :: Problem.AbstractProblem,
                         time :: Float64)
    
    n = ref_el.n_node
    d_loc = zeros(n, n, mesh.n_cell)
    
    
    x = FEM.map_ref_point(dof, quad.point, 1:dof.n_elem)
    diffusion = Problem.diffusion(problem, time, x)

    DF = FEM.map_ref_point_grad_inv(dof, quad.point, 1:dof.n_elem);
    weight_elem = FEM.map_ref_point_grad_det(dof, quad.point, 1:dof.n_elem)

    DPhi = FEM.shapeFun_grad(ref_el, quad.point)
    DPhi_test = FEM.shapeFun_grad(ref_el, quad.point)

    for k = 1:mesh.n_cell
    	for j=1:length(quad.weight)
    		for i=1:ref_el.n_node
        		d_loc[i,j,k] = weight_elem[k] * quad.weight[j] * (DPhi_test[i,j]'*DF[k]) * (diffusion[k][j]*(DF[k]'*DPhi[i,j]))
        	end
        end
    end

    return -d_loc
end # end function


# -----------------------------
# -----------------------------