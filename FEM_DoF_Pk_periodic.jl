type Dof_Pk_periodic_square <: Dof_square

    mesh :: Mesh.TriMesh
    
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
    # Topology info and related functions
    is_periodic :: Bool
    n_true_dof :: Int64

    map_vec_ind_mesh2dof :: Array{Int64,1}
    #map_vec_ind_dof2mesh :: Array{Int64,1}

    map_ind_dof2mesh :: Function
    # ----------------------------------------

    get_dof_elem :: Function
    get_dof_edge :: Function

    function Dof_Pk_periodic_square(mesh :: Mesh.TriMesh,
                             ref_el :: RefEl_Pk)
        
        this = new()
            
        # ----------------------------------------
        this.FEM_order = ref_el.n_order
        this.FEM_info = "DoF object for ---   periodic unit square  --- Pk-Lagrange FEM of order $(ref_el.n_order)."
        # ----------------------------------------

        if ref_el.n_order==1
            # ----------------------------------------------------------------------------------------------------------------------------------------

            this.mesh = mesh
            
            # ----------------------------------------
            # Node infos
            this.n_node = mesh.n_point
            this.n_node_boundary = 0
            this.n_node_interior = mesh.n_point
            this.n_node_dirichlet = 0
            this.n_node_neumann = 0
            this.n_node_per_edge = 2
            this.n_node_per_elem = 3
            
            this.ind_node_boundary = []
            this.ind_node_interior = collect(1:mesh.n_point)
            this.ind_node_dirichlet  = []
            this.ind_node_non_dirichlet = setdiff(1:this.n_node, this.ind_node_dirichlet)
            this.ind_node_neumann = []
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Edge infos
            this.n_edge = mesh.n_edge
            this.n_edge_boundary = 0
            this.n_edge_interior = mesh.n_edge
            #this.n_edge_dirichlet = []
            #this.n_edge_neumann = []
            this.n_edge_per_elem = 3
            
            this.ind_edge_boundary  = []
            this.ind_edge_interior = collect(1:mesh.n_edge)
            #this.ind_edge_dirichlet = []
            #this.ind_edge_neumann = []
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Element infos
            this.n_elem = mesh.n_cell
            this.n_elem_boundary = 0
            this.n_elem_interior = mesh.n_cell
            
            this.ind_elem_boundary = []
            this.ind_elem_interior = collect(1:mesh.n_cell)
            # ----------------------------------------
            
            
            # ----------------------------------------
            # Topology info
            this.is_periodic = true

            n_node_boundary = sum(mesh.point_marker.==1)
            this.n_true_dof = mesh.n_point - 3 - 0.5*(n_node_boundary - 4)
            # ----------------------------------------

            
            # ----------------------------------------
            # Create the maps that translate dofs to mesh based
            # variables. This is one of the core modules for the
            # handling DoFs of periodic meshes (without actually
            # touching the mesh).

            # Tells where any mesh based variable can be found
            # topologically
            this.map_vec_ind_mesh2dof = zeros(mesh.n_point)

            n_edge_per_seg = Int64(size(mesh.segment,1)/4)
            
            # Corners of the square
            this.map_vec_ind_mesh2dof[1] = 1
            this.map_vec_ind_mesh2dof[1+n_edge_per_seg] = 1
            this.map_vec_ind_mesh2dof[1+2*n_edge_per_seg] = 1
            this.map_vec_ind_mesh2dof[1+3*n_edge_per_seg] = 1

            # Segment 1 (0,0)->(1,0)
            for i=1:(n_edge_per_seg-1)
                this.map_vec_ind_mesh2dof[1+i] = 1 + i
            end

            # Segment 2 (1,0)->(1,1)
            for i=1:(n_edge_per_seg-1)
                this.map_vec_ind_mesh2dof[1+n_edge_per_seg+i] = 1 + n_edge_per_seg + i - 1
            end

            # Segment 3 (1,1)->(0,1)
            for i=1:(n_edge_per_seg-1)
                this.map_vec_ind_mesh2dof[1+2*n_edge_per_seg+i] = this.map_vec_ind_mesh2dof[1+n_edge_per_seg - i]
            end

            # Segment 4 (0,1)->(0,0)
            for i=1:(n_edge_per_seg-1)
                this.map_vec_ind_mesh2dof[1+3*n_edge_per_seg+i] = this.map_vec_ind_mesh2dof[1+2*n_edge_per_seg - i]
            end

            this.map_vec_ind_mesh2dof[n_node_boundary+1:end] = collect(   ((n_node_boundary + 1):(mesh.n_point)) - 3 - 0.5*(n_node_boundary - 4)  )
            
            this.map_ind_dof2mesh = function(vec_dof :: Array{Float64})
                # Map a vector in terms of degrees of freedom to a
                # vector on the actual mesh. This is only interesting
                # for back mapping. The map mesh->dof can be found in
                # the function get_dof()
                error("Check the mappings!!!")
                vec_mesh = vec_dof[this.map_vec_ind_mesh2dof[sortperm(this.map_vec_ind_mesh2dof)]][sortperm(sortperm(this.map_vec_ind_mesh2dof))]
            end
            # ----------------------------------------
            
            
            # ----------------------------------------
            this.get_dof_elem = function(ind_c)
                ind_c = vec(collect(ind_c))

                # Put a filter before mesh indices to translate to dof
                # indices
                dofs = this.map_vec_ind_mesh2dof[this.mesh.get_cell(ind_c)]
            
                return dofs
            end # end function


            # ----------------------------------------
            this.get_dof_edge = function(ind_e)
                ind_c = vec(collect(ind_e))

                # Put a filter before mesh indices to translate to dof
                # indices
                dofs = this.map_vec_ind_mesh2dof[this.mesh.get_edge(ind_c)]
            
                return dofs
            end # end function
            # ----------------------------------------
            # ----------------------------------------------------------------------------------------------------------------------------------------
                
        elseif ref_el.n_order==2
            error("Order 2 DoF not implemented yet (periodic setting).")
        end # end if order
        
        return this    
    end # end constructor
end # end type
