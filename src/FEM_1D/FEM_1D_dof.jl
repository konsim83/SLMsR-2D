struct Dof_1D
    order :: Real

    n_mesh_elem :: Int64
    n_mesh_point :: Int64

    n_nodes :: Int64
    n_nodes_boundary :: Int64
    n_nodes_interior :: Int64
    n_nodes_dirichlet :: Int64

    n_nodes_per_elem :: Int64

    n_true_dof :: Int64

    ind_trial :: Array{Int,1}
    ind_test :: Array{Int,1}

    function Dof_1D(mesh :: Mesh_1D)
        
        order = 1

        n_mesh_point = mesh.n_point
        n_mesh_elem = mesh.n_cell

        n_nodes = n_mesh_point
        n_nodes_boundary = 2
        n_nodes_interior = n_nodes - 2
        n_nodes_dirichlet = 2
        
        n_nodes_per_elem = 2
        
        n_true_dof = n_nodes

        ind_test = vec([mesh.cell ; mesh.cell])
        ind_trial = vec(transpose([mesh.cell[1,:] mesh.cell[1,:] mesh.cell[2,:] mesh.cell[2,:]]))
        
        return new(order,
                    n_mesh_elem,
                    n_mesh_point,
                    n_nodes,
                    n_nodes_boundary,
                    n_nodes_interior,
                    n_nodes_dirichlet,
                    n_nodes_per_elem,
                    n_true_dof,
                    ind_trial, 
                    ind_test)
    end # end constructor
end # end type

function get_dofs(dof :: Dof_1D, mesh :: Mesh_1D, ind_c :: Array{Int64,1})
    
    return mesh.cell[:,ind_c]
end # end function