type Dof_Pk_periodic
    
    # ----------------------------------------
    # General infos
    FEM_order :: Int64
    FEM_info :: Int64
    # ----------------------------------------


    # ----------------------------------------
    # Node infos
    n_node :: Int64
    n_node_boundary :: Int64
    n_node_interior :: Int64
    #n_node_dirichlet :: Int64
    #n_node_neumann :: Int64
    n_node_per_edge :: Int64
    n_node_per_elem :: Int64

    ind_node_boundary :: Array{Int64,1}
    ind_node_interior :: Array{Int64,1}
    #ind_node_dirichlet :: Array{Int64,1}
    #ind_node_neumann :: Array{Int64,1}
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
    # If the topology is periodic
    n_true_dof :: Int64
    # ----------------------------------------


    get_dofs :: Function

    function Dof_Pk_periodic(mesh :: TriMesh,
                             ref_el :: RefEl_Pk)
        
        this = new()
            
        # ----------------------------------------
        this.FEM_order = ref_el.order
        this.FEM_order = "DoF object for ---   periodic   --- Pk-Lagrange FEM of order $ref_el.order."
        # ----------------------------------------

        if ref_el.order==1
            # ----------------------------------------------------------------------------------------------------------------------------------------
            # ----------------------------------------
            # Node infos
            this.n_node = 
                this.n_node_boundary :: Int64
            this.n_node_interior :: Int64
            #this.n_node_dirichlet = []
            #this.n_node_neumann = []
            this.n_node_per_edge :: Int64
            this.n_node_per_elem :: Int64
            
            this.ind_node_boundary :: Array{Int64,1}
            this.ind_node_interior :: Array{Int64,1}
            #this.ind_node_dirichlet  = []
            #this.ind_node_neumann = []
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Edge infos
            this.n_edge :: Int64
            this.n_edge_boundary :: Int64
            this.n_edge_interior :: Int64
            #this.n_edge_dirichlet = []
            #this.n_edge_neumann = []
            this.n_edge_per_elem :: Int64
            
            this.ind_edge_boundary :: Int64
            this.ind_edge_interior :: Int64
            #this.ind_edge_dirichlet = []
            #this.ind_edge_neumann = []
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Element infos
            this.n_elem :: Int64
            this.n_elem_boundary :: Int64
            this.n_elem_interior :: Int64
            
            this.ind_elem_boundary :: Array{Int64,1}
            this.ind_elem_interior :: Array{Int64,1}
            # ----------------------------------------
            
            
            # ----------------------------------------
            this.get_dofs = function(ind_c)
                
                ind_c = vec(collect(ind_c))
                
                ind_c = (this.order) * repmat(ind_c - 1, 1, this.order+1) + 1
                ind_loc = transpose(   repmat(collect(1:this.order+1) , 1, size(ind_c,1))   ) - 1
                d = ind_c + ind_loc
                if this.is_periodic
                    d[collect(1:length(d))[collect(d.==this.n_nodes_all)]] = 1
                end
            
                return d
            end # end function
            # ----------------------------------------
            # ----------------------------------------------------------------------------------------------------------------------------------------
                
        elseif ref_el.order==2
            error("Order 2 DoF not implemented yet (periodic setting).")
        end # end if order
        
        return this    
    end # end constructor
end # end type
