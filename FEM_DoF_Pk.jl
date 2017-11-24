type Dof_Pk{FEM_order} <: Dof
    
    # ----------------------------------------
    # General infos
    FEM_order :: Int64
    FEM_info :: String
    # ----------------------------------------


    # ----------------------------------------
    # Node infos
    n_node :: Int64
    n_node_boundary :: Int64
    n_node_interior :: Int64
    n_node_dirichlet :: Int64
    n_node_neumann :: Int64
    n_node_per_edge :: Int64
    n_node_per_elem :: Int64

    ind_node_boundary :: Array{Int64,1}
    ind_node_interior :: Array{Int64,1}
    ind_node_dirichlet :: Array{Int64,1}
    ind_node_non_dirichlet :: Array{Int64,1}
    ind_node_neumann :: Array{Int64,1}
    # ----------------------------------------


    # ----------------------------------------
    # Edge infos
    n_edge :: Int64
    n_edge_boundary :: Int64
    n_edge_interior :: Int64
    #n_edge_dirichlet :: Int64
    #n_edge_neumann :: Int64
    n_edge_per_elem :: Int64

    ind_edge_boundary :: Array{Int64,1}
    ind_edge_interior :: Array{Int64,1}
    #ind_edge_dirichlet :: Array{Int64,1}
    #ind_edge_neumann :: Array{Int64,1}
    # ----------------------------------------

    # ----------------------------------------
    # Element infos
    n_elem :: Int64
    n_elem_boundary :: Int64
    n_elem_interior :: Int64

    ind_elem_boundary :: Array{Int64,1}
    ind_elem_interior :: Array{Int64,1}
    # ----------------------------------------
    
    # ----------------------------------------
    # Topology info
    is_periodic :: Bool
    n_true_dof :: Int64

    map_vec_ind_mesh2dof :: Array{Int64,1}
    #map_vec_ind_dof2mesh :: Array{Int64,1}
    # ----------------------------------------
    

    function Dof_Pk(mesh :: Mesh.TriMesh)
        
        this = new()
            
        # ----------------------------------------
        this.FEM_order = FEM_order 
        this.FEM_info = "DoF object for ---   non-periodic   --- Pk-Lagrange FEM of order $(FEM_order)."
        # ----------------------------------------

        if FEM_order==1
            # ----------------------------------------------------------------------------------------------------------------------------------------
            # ----------------------------------------
            # Node infos
            this.n_node = mesh.n_point
            this.n_node_boundary = sum(mesh.point_marker.==1)
            this.n_node_interior = sum(mesh.point_marker.==0)
            this.n_node_dirichlet = this.n_node_boundary
            this.n_node_neumann = 0
            this.n_node_per_edge = 2
            this.n_node_per_elem = 3
            
            this.ind_node_boundary = find(mesh.point_marker.==1)
            this.ind_node_interior = find(mesh.point_marker.==0)
            this.ind_node_dirichlet  = find(mesh.point_marker.==1)
            this.ind_node_non_dirichlet = setdiff(1:this.n_node, this.ind_node_dirichlet)
            this.ind_node_neumann  = []
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Edge infos
            this.n_edge = mesh.n_edge
            this.n_edge_boundary = sum(mesh.edge_marker.==1)
            this.n_edge_interior = sum(mesh.edge_marker.==0)
            #this.n_edge_dirichlet = []
            #this.n_edge_neumann = []
            this.n_edge_per_elem = 3
            
            this.ind_edge_boundary = find(mesh.edge_marker.==1)
            this.ind_edge_interior = find(mesh.edge_marker.==0)
            #this.ind_edge_dirichlet = []
            #this.ind_edge_neumann = []
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Element infos
            this.n_elem = mesh.n_cell
            this.n_elem_boundary = sum(sum(mesh.cell_neighbor.==0,2).!=0)
            this.n_elem_interior = sum(sum(mesh.cell_neighbor.==0,2).==0)
            
            this.ind_elem_boundary = find(sum(mesh.cell_neighbor.==0,2).!=0)
            this.ind_elem_interior = find(sum(mesh.cell_neighbor.==0,2).==0)
            # ----------------------------------------

            
            # ----------------------------------------
            # Topology info
            this.is_periodic = false
            this.n_true_dof = mesh.n_point

            this.map_vec_ind_mesh2dof = collect(1:this.n_node)
            #this.map_vec_ind_dof2mesh :: Array{Int64,1}
            # ----------------------------------------
            # ----------------------------------------------------------------------------------------------------------------------------------------

        elseif FEM_order==2
            error("Order 2 DoF not implemented yet.")
        end # end if order
        
        return this
    end # end constructor
end # end type



# -------------------------------------------------------------------------------------------
# -----------------------------   Functions on Dof_Pk   -----------------------------
# -------------------------------------------------------------------------------------------

# ----------------------------------------
function map_ind_dof2mesh(dof :: Dof_Pk{1}, vec_dof :: Array{Float64})
    """

    Map a vector in terms of degrees of freedom to a vector on the
    actual mesh. This is only interesting for back mapping. Here
    the mapping is trivial

    """

    return vec_dof
end
# ----------------------------------------


# ----------------------------------------
function get_dof_elem(dof :: Dof_Pk{1}, mesh :: Mesh.TriMesh, ind_c :: Array{Int64,1})
    
    return Mesh.get_cell(mesh, ind_c)
end # end function

function get_dof_elem(dof :: Dof_Pk{1}, mesh :: Mesh.TriMesh, ind_c :: UnitRange{Int64})
    
    return Mesh.get_cell(mesh, ind_c)
end # end function
# ----------------------------------------


# ----------------------------------------
function get_dof_edge(dof :: Dof_Pk{1}, mesh :: Mesh.TriMesh, ind_e :: Array{Int64,1})
    
    return Mesh.get_edge(mesh, ind_e)
end # end function

function get_dof_edge(dof :: Dof_Pk{1}, mesh :: Mesh.TriMesh, ind_e :: UnitRange{Int64})
    
    return Mesh.get_edge(mesh, ind_e)
end # end function
# ----------------------------------------

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------
