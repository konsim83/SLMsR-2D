type Dof
    order :: Real

    n_mesh_elem :: Int64
    n_mesh_point :: Int64

    is_periodic :: Bool

    n_nodes_all :: Int64
    n_nodes_boundary :: Int64
    n_true_dof :: Int64
    n_nodes_interior :: Int64
    n_nodes_dirichlet :: Int64
    n_nodes_neumann :: Int64
    n_bubble_node :: Int64
    n_nodes_per_elem :: Int64
    n_quad_points_all :: Int64

    nodes_all
    quad_points_all

    get_dofs :: Function

    function Dof(mesh::Mesh,
                 ref_el::RefEl_Lagrange,
                 quad::Quad_jacobi_gauss,
                 is_periodic)
        
        this = new()
        
        this.order = ref_el.order

        this.n_mesh_point = mesh.n_point
        this.n_mesh_elem = mesh.n_elem

        this.n_nodes_all = this.n_mesh_elem*(this.order - 1) + this.n_mesh_point
        this.is_periodic = is_periodic
        this.n_true_dof = this.n_nodes_all - this.is_periodic
        this.n_nodes_boundary = 2*this.is_periodic            
        this.n_nodes_interior = this.n_nodes_all - this.n_nodes_boundary
        this.n_nodes_dirichlet = -1
        this.n_nodes_neumann = -1
        this.n_bubble_node = this.n_mesh_elem * (this.order - 1)
        this.n_nodes_per_elem = this.order + 1

        # This saves operations
        this.nodes_all = push!(vec(transpose(map_ref_point(mesh, ref_el.node[1:end-1]))), mesh.b)
        this.quad_points_all = vec(map_ref_point(mesh, quad.point)')
        this.n_quad_points_all = length(this.quad_points_all)
        

        #this.ind_nodes_boundary = []
        #this.ind_nodes_interior = []
        #this.ind_nodes_dirichlet = []
        #this.ind_nodes_neumann = []
        #this.ind_bubble_node = []


        this.get_dofs = function(ind_c)
            
            ind_c = vec(collect(ind_c))
            
            ind_c = (this.order) * repmat(ind_c - 1, 1, this.order+1) + 1
            ind_loc = transpose(   repmat(collect(1:this.order+1) , 1, size(ind_c,1))   ) - 1
            d = ind_c + ind_loc
            if this.is_periodic
                d[find(d.==this.n_nodes_all)] = 1
            end
            
            return d
        end # end function
        
        
        return this
    end # end constructor
end # end type


type Dof_reconstruction
    order :: Real

    n_mesh_elem :: Int64
    n_mesh_point :: Int64

    is_periodic :: Bool

    n_nodes_all :: Int64
    n_nodes_boundary :: Int64
    n_true_dof :: Int64
    n_nodes_interior :: Int64
    n_nodes_dirichlet :: Int64
    n_nodes_neumann :: Int64
    n_bubble_node :: Int64
    n_nodes_per_elem :: Int64
    n_quad_points_all :: Int64

    nodes_all
    quad_points_all

    #ind_nodes_boundary
    #ind_nodes_interior
    #ind_nodes_dirichlet
    #ind_nodes_neumann
    #ind_bubble_node

    get_dofs :: Function
    get_dofs_eval :: Function

    function Dof_reconstruction(mesh::Mesh,
                 ref_el::RefEl_Lagrange,
                 quad::Quad_jacobi_gauss,
                 is_periodic)
        
        this = new()
        
        this.order = ref_el.order

        this.n_mesh_point = mesh.n_point
        this.n_mesh_elem = mesh.n_elem

        this.n_nodes_all = this.n_mesh_elem*(this.order - 1) + this.n_mesh_point
        this.is_periodic = is_periodic
        this.n_true_dof = (this.n_nodes_all - this.is_periodic) + mesh.n_elem*2
        this.n_nodes_boundary = 2*this.is_periodic            
        this.n_nodes_interior = this.n_nodes_all - this.n_nodes_boundary
        this.n_nodes_dirichlet = -1
        this.n_nodes_neumann = -1
        this.n_bubble_node = this.n_mesh_elem * (this.order - 1)
        this.n_nodes_per_elem = this.order + 1

        # This saves operations
        this.nodes_all = push!(vec(transpose(map_ref_point(mesh, ref_el.node[1:end-1]))), mesh.b)
        this.quad_points_all = vec(map_ref_point(mesh, quad.point)')
        this.n_quad_points_all = length(this.quad_points_all)
        

        #this.ind_nodes_boundary = []
        #this.ind_nodes_interior = []
        #this.ind_nodes_dirichlet = []
        #this.ind_nodes_neumann = []
        #this.ind_bubble_node = []


        this.get_dofs = function(ind_c)
            
            ind_c = vec(collect(ind_c))
            ind_c_orig = vec(collect(ind_c))
            
            ind_c = (this.order) * repmat(ind_c - 1, 1, this.order+1) + 1
            ind_loc = transpose(   repmat(collect(1:this.order+1) , 1, size(ind_c,1))   ) - 1
            d = ind_c + ind_loc
            
            if this.is_periodic
                d[find(d.==this.n_nodes_all)] = 1
            end
            
            # display([d d+this.n_nodes_all-this.is_periodic])

            return hcat(d, [2*ind_c_orig-1 2*ind_c_orig]+this.n_nodes_all-this.is_periodic)
        end # end function


        this.get_dofs_eval = function(ind_c)
            
            ind_c = vec(collect(ind_c))
            ind_c_orig = vec(collect(ind_c))
            
            ind_c = (this.order) * repmat(ind_c - 1, 1, this.order+1) + 1
            ind_loc = transpose(   repmat(collect(1:this.order+1) , 1, size(ind_c,1))   ) - 1
            d = ind_c + ind_loc
            
            if this.is_periodic
                d[find(d.==this.n_nodes_all)] = 1
            end
            
            # display([d d+this.n_nodes_all-this.is_periodic])

            return d, [2*ind_c_orig-1 2*ind_c_orig]+this.n_nodes_all-this.is_periodic #d+this.n_nodes_all-this.is_periodic
        end # end function
        
        
        return this
    end # end constructor
end # end type