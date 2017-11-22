type Dof_Pk <: Dof
    
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
    # ----------------------------------------


    get_dof_elem :: Function
    get_dof_edge :: Function

    function Dof_Pk(mesh :: Mesh.TriMesh,
                    ref_el :: RefEl_Pk)
        
        this = new()
            
        # ----------------------------------------
        this.FEM_order = ref_el.n_order
        this.FEM_info = "DoF object for ---   non-periodic   --- Pk-Lagrange FEM of order $(ref_el.n_order)."
        # ----------------------------------------

        if ref_el.n_order==1
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
            # ----------------------------------------
            
            
            # ----------------------------------------
            # P1 version
            this.get_dof_elem = function(ind_c)
                dofs = mesh.get_cell(ind_c)
                
                return dofs
            end # end function


            this.get_dof_edge = function(ind_e)
                dofs = mesh.get_edge(ind_e)
                
                return dofs
            end # end function
            # ----------------------------------------
            # ----------------------------------------------------------------------------------------------------------------------------------------

        elseif ref_el.n_order==2
            error("Order 2 DoF not implemented yet.")
        end # end if order
        
        return this
    end # end constructor
end # end type
