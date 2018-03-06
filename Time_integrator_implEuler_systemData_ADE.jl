# ------------------------------------------------------------------------------
# ------- System data   -------
struct System_data_implEuler_ADE <: AbstractSystem_data_implEuler

    mass :: SparseMatrixCSC{Float64,Int64}
    advection :: SparseMatrixCSC{Float64,Int64}
    diffusion :: SparseMatrixCSC{Float64,Int64}
    rhs :: Array{Float64}
    
    system_matrix :: SparseMatrixCSC{Float64,Int64}
    system_rhs :: Array{Float64}

    # holds solution of system in dof-indices. Note that Dirichlet
    # nodes can occur with periodic boundaries.
    u_temp :: Array{Float64}
        
end


function System_data_implEuler_ADE(dof :: FEM.AbstractDof, mesh :: Mesh.TriangleMesh.TriMesh, problem :: Problem.AbstractPhysicalProblem)
        
    # Create a pattern
    i = FEM.get_dof_elem(dof, mesh, 1:dof.n_elem)
    ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
    ind_test = vec(transpose(repmat(i, 1, size(i,2))))

    # Allocate memory for sparse matrix
    mass = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    advection = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    diffusion = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    rhs = zeros(dof.n_true_dof)

    system_matrix = mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
    system_rhs = rhs[dof.ind_node_non_dirichlet]

    # initialize a temporary array that holds the 
    u_temp = FEM.map_vec_mesh2dof(dof, Problem.u_init(problem, mesh.point))
   
    return System_data_implEuler_ADE(mass, advection, diffusion, rhs, system_matrix, system_rhs, u_temp)
end


function System_data_implEuler_ADE(dof :: FEM.AbstractDof, mesh :: Mesh.TriangleMesh.TriMesh, problem :: Problem.AbstractBasisProblem)
    
    # Create a pattern
    i = FEM.get_dof_elem(dof, mesh, 1:dof.n_elem)
    ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
    ind_test = vec(transpose(repmat(i, 1, size(i,2))))

    # Allocate memory for sparse matrix
    mass = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    advection = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    diffusion = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
    rhs = zeros(dof.n_true_dof,3)

    system_matrix = mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
    system_rhs = rhs[dof.ind_node_non_dirichlet,:]

    # Initialize a temporary array that holds the previous time
    # step data. Dirichlet data if fixed over iterations. Used for
    # update.
    u_temp = FEM.map_vec_mesh2dof(dof, Problem.u_init(problem, mesh.point))
   
    return System_data_implEuler_ADE(mass, advection, diffusion, rhs, system_matrix, system_rhs, u_temp)
end
# ------------------------------------------------------------------------------
