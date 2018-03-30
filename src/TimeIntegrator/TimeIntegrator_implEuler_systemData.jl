# ------------------------------------------------------------------------------
# ------- System data   -------
struct ImplEulerData <: AbstractSystemData

    # Solves M*u_New = M*uOld + A*uNew + rhs

    M :: SparseMatrixCSC{Float64,Int64}
    A :: SparseMatrixCSC{Float64,Int64}
    rhs :: Array{Float64}
    
    system_matrix :: SparseMatrixCSC{Float64,Int64}
    system_rhs :: Array{Float64}

    # Contais topological indices of Dirichlet and non-Dirichlet nodes.
    ind_node_non_dirichlet :: Array{Int,1}
    ind_node_dirichlet :: Array{Int,1}
        
end


# Constructor for physical problem
function ImplEulerData(dof :: FEM.AbstractDof,
                        mesh :: Mesh.TriangleMesh.TriMesh,
                        problem :: Problem.AbstractPhysicalProblem)

    # ------------------------------
    # Map pure mesh indices to dof indices
    ind_node_non_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_non_dirichlet))
    ind_node_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_dirichlet))
    # ------------------------------

    # ------------------------------
    # Allocate memory for sparse matrix        
    ind = dof.ind
    ind_test = dof.ind_test

    M = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    A = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    rhs = zeros(dof.n_true_dof)
    # ------------------------------

    system_matrix = M[ind_node_non_dirichlet,ind_node_non_dirichlet]
    system_rhs = rhs[ind_node_non_dirichlet]
   
    return ImplEulerData(M, A, rhs,
                            system_matrix, system_rhs, 
                            ind_node_non_dirichlet, ind_node_dirichlet)
end


# Constructor for basis problem
function ImplEulerData(dof :: FEM.AbstractDof, 
                            mesh :: Mesh.TriangleMesh.TriMesh, 
                            problem :: Problem.AbstractBasisProblem)
    
    # ------------------------------
    ind_node_non_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_non_dirichlet))
    ind_node_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_dirichlet))
    # ------------------------------

    # ------------------------------
    # Allocate memory for sparse matrix        
    ind = dof.ind
    ind_test = dof.ind_test

    M = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    A = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    rhs = zeros(dof.n_true_dof,3)
    # ------------------------------

    system_matrix = mass[ind_node_non_dirichlet,ind_node_non_dirichlet]
    system_rhs = rhs[ind_node_non_dirichlet,:]
   
    return ImplEulerData(M, A, rhs,
                            system_matrix, system_rhs, 
                            ind_node_non_dirichlet, ind_node_dirichlet)
end
# ------------------------------------------------------------------------------
