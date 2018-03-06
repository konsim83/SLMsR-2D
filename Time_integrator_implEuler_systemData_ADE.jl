# ------------------------------------------------------------------------------
# ------- System data   -------
struct System_data_implEuler_ADE <: AbstractSystem_data_implEuler

    mass :: SparseMatrixCSC{Float64,Int64}
    advection :: SparseMatrixCSC{Float64,Int64}
    diffusion :: SparseMatrixCSC{Float64,Int64}
    rhs :: Array{Float64}
    
    system_matrix :: SparseMatrixCSC{Float64,Int64}
    system_rhs :: Array{Float64}

    # Holds solution of system in dof-indices. Note that Dirichlet
    # nodes can occur with periodic boundaries.
    u_temp :: Array{Float64}

    # Contais topological indices of Dirichlet and non-Dirichlet nodes.
    ind_node_non_dirichlet :: Array{Int,1}
    ind_node_dirichlet :: Array{Int,1}
        
end


function System_data_implEuler_ADE(dof :: FEM.AbstractDof, mesh :: Mesh.TriangleMesh.TriMesh, problem :: Problem.AbstractPhysicalProblem)

    # ------------------------------
    ind_node_non_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_non_dirichlet))
    ind_node_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_dirichlet))
    # ------------------------------

    # ------------------------------
    # Allocate memory for sparse matrix        
    ind = dof.ind
    ind_test = dof.ind_test

    mass = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    advection = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    diffusion = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    rhs = zeros(dof.n_true_dof)
    # ------------------------------

    system_matrix = mass[ind_node_non_dirichlet,ind_node_non_dirichlet]
    system_rhs = rhs[ind_node_non_dirichlet]

    # initialize a temporary array that holds the data
    u_temp = FEM.map_vec_mesh2dof(dof, Problem.u_init(problem, mesh.point))
   
    return System_data_implEuler_ADE(mass, advection, diffusion, rhs,
                                        system_matrix, system_rhs, 
                                        u_temp, 
                                        ind_node_non_dirichlet, ind_node_dirichlet)
end


function System_data_implEuler_ADE(dof :: FEM.AbstractDof, mesh :: Mesh.TriangleMesh.TriMesh, problem :: Problem.AbstractBasisProblem)
    
    # ------------------------------
    ind_node_non_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_non_dirichlet))
    ind_node_dirichlet = unique(FEM.map_ind_mesh2dof(dof, dof.ind_node_dirichlet))
    # ------------------------------

    # ------------------------------
    # Allocate memory for sparse matrix        
    ind = dof.ind
    ind_test = dof.ind_test

    mass = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    advection = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    diffusion = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    rhs = zeros(dof.n_true_dof,3)
    # ------------------------------

    system_matrix = mass[ind_node_non_dirichlet,ind_node_non_dirichlet]
    system_rhs = rhs[ind_node_non_dirichlet,:]

    # Initialize a temporary array that holds the previous time
    # step data. Dirichlet data if fixed over iterations. Used for
    # update.
    u_temp = FEM.map_vec_mesh2dof(dof, Problem.u_init(problem, mesh.point))
   
    return System_data_implEuler_ADE(mass, advection, diffusion, rhs,
                                        system_matrix, system_rhs, 
                                        u_temp, 
                                        ind_node_non_dirichlet, ind_node_dirichlet)
end
# ------------------------------------------------------------------------------
