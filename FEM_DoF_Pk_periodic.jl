type Dof_Pk_periodic_square{FEM_order} <: AbstractDof
    
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

    index_map_dof2mesh :: Array{Int64,1}
    index_map_mesh2dof :: Array{Int64,1}
    # ----------------------------------------

    ind :: Array{Int64,1}
    ind_test :: Array{Int64,1}
    ind_lin :: Array{Int64,1}    
    
    T_ref2cell :: Array{Float64,3}
    T_cell2ref :: Array{Float64,3}

    function Dof_Pk_periodic_square{FEM_order}(mesh :: Mesh.TriangleMesh.TriMesh) where {FEM_order}
        
        this = new()
            
        # ----------------------------------------
        this.FEM_order = FEM_order
        this.FEM_info = "DoF object for ---   periodic unit square  --- Pk-Lagrange FEM of order $(FEM_order)."
        # ----------------------------------------

        if FEM_order==1
            # ----------------------------------------------------------------------------------------------------------------------------------------
            

            # ----------------------------------------
            # Create the maps that translate dofs to mesh based
            # variables. This is one of the core modules for the
            # handling DoFs of periodic meshes (without actually
            # touching the mesh).

            # Tells where any mesh based variable can be found
            # topologically

            # Needs information on  the mesh
            edge_marker_pair = [1 3 ; 2 4]
            f_edge_to_edge = Array{Function,1}(2)
            f_edge_to_edge[1] = function(p :: Array{Float64,2})
                # Periodicity in y-direction, maps upper edge to lower edge

                return broadcast(+, [0.0 -1.0], p)
            end

            f_edge_to_edge[2] = function(p :: Array{Float64,2})
                # Periodicity in x-direction, maps left edge to right edge

                return broadcast(+, [1.0 0.0], p)
            end
            this.index_map_dof2mesh = identify_points(mesh, edge_marker_pair, f_edge_to_edge)
            this.index_map_mesh2dof = indexin(unique(this.index_map_dof2mesh), this.index_map_dof2mesh)

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

            n_node_boundary = sum(mesh.point_marker.!=0)
            n_true_dof = Int(mesh.n_point - 3 - 0.5*(n_node_boundary - 4))
            this.ind_node_non_dirichlet = setdiff(1:n_true_dof, this.ind_node_dirichlet)
            
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

            this.n_true_dof = n_true_dof
            # ----------------------------------------

            ind_cell = this.index_map_dof2mesh[Mesh.get_cell(mesh, 1:mesh.n_cell)]   
            this.ind = vec(ind_cell[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
            this.ind_test = vec(transpose(repmat(ind_cell, 1, size(ind_cell,2))))
            this.ind_lin = sub2ind((this.n_true_dof,this.n_true_dof), this.ind_test, this.ind)
            
            this.T_ref2cell = zeros(2, 3, mesh.n_cell)
            this.T_cell2ref = zeros(2, 3, mesh.n_cell)
            for i=1:this.n_elem
                # Extended transformation matrices for mapping from and to the
                # reference cell K = [(0,0), (1,0), (0,1)]
                this.T_ref2cell[:,:,i] = [   mesh.point[mesh.cell[i,2],:]-mesh.point[mesh.cell[i,1],:]   mesh.point[mesh.cell[i,3],:]-mesh.point[mesh.cell[i,1],:]   mesh.point[mesh.cell[i,1],:]   ]
                this.T_cell2ref[:,:,i] = (eye(2)/this.T_ref2cell[:,1:2,i]) * [   eye(2)  mesh.point[mesh.cell[i,1],:]   ]
            end

            # ----------------------------------------------------------------------------------------------------------------------------------------
                
        elseif FEM_order==2
            error("Order 2 DoF not implemented yet (periodic setting).")
        end # end if order
        
        return this    
    end # end constructor
end # end type



# -------------------------------------------------------------------------------------------------------
# -----------------------------   Functions on Dof_Pk_periodic   -----------------------------
# -------------------------------------------------------------------------------------------------------

# ----------------------------------------
function map_ind_dof2mesh(dof :: Dof_Pk_periodic_square{1}, vec_dof :: Array{Float64})

    """

    Map a vector in terms of degrees of freedom to a vector on the
    actual mesh. This is only interesting for back mapping. The map
    mesh->dof can be found in the function get_dof()

    """

    # expand a dof-vector into a mesh-vector (only periodic boundaries)
    vec_mesh = vec_dof[dof.index_map_dof2mesh,:]

    return vec_mesh
end
# ----------------------------------------


# ----------------------------------------
function map_ind_mesh2dof(dof :: Dof_Pk_periodic_square{1}, vec_mesh :: Array{Float64})

    """

    Map a vector of Float values in terms of node indices to degrees of
    freedom indices.

    """

    # reduce a mesh-vector to a dof-vector (only periodic boundaries)
    vec_dof = vec_mesh[dof.index_map_mesh2dof,:]

    return vec_dof
end
# ----------------------------------------


# ----------------------------------------
function get_dof_elem(dof :: Dof_Pk_periodic_square{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_c :: Array{Int64,1})
    # Put a filter before mesh indices to translate to dof
    # indices
    dofs = dof.index_map_dof2mesh[Mesh.get_cell(mesh, ind_c)]
    
    return dofs
end # end function

function get_dof_elem(dof :: Dof_Pk_periodic_square{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_c :: UnitRange{Int64})
    ind_c = vec(collect(ind_c))
    
    # Put a filter before mesh indices to translate to dof
    # indices
    dofs = dof.index_map_dof2mesh[Mesh.get_cell(mesh, ind_c)]
    
    return dofs
end # end function
# ----------------------------------------


# ----------------------------------------
function get_dof_edge(dof :: Dof_Pk_periodic_square{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_e :: Array{Int64,1})    
    # Put a filter before mesh indices to translate to dof
    # indices
    dofs = dof.index_map_dof2mesh[Mesh.get_edge(mesh, ind_c)]
    
    return dofs
end # end function

function get_dof_edge(dof :: Dof_Pk_periodic_square{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_e :: UnitRange{Int64})
    ind_c = vec(collect(ind_e))
    
    # Put a filter before mesh indices to translate to dof
    # indices
    dofs = dof.index_map_dof2mesh[Mesh.get_edge(mesh, ind_c)]
    
    return dofs
end # end function
# ----------------------------------------


# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------


# ----------------------------------------
function identify_points(mesh :: Mesh.TriangleMesh.TriMesh,
                            edge_marker_pair :: Array{Int64,2},
                            edge_trafo :: Array{Function,1})

    length(edge_trafo)!=size(edge_marker_pair,1) ? 
        error("Number of Matched edge pairs must coincide with transformations.") :


    ind_point_boundary = find(mesh.point_marker.!=0)
    ind_point_inner = find(mesh.point_marker.==0)

    index_map = Array{Array{Int64,1},1}(0)

    for i in ind_point_boundary
        push!(index_map,[i])
    end

    for i in 1:size(edge_marker_pair,1)
        pair = edge_marker_pair[i,:]

        point_ind_on_edge1 = sort(unique(vec(mesh.edge[mesh.edge_marker.==pair[1],:]')))
        # display(point_ind_on_edge1)
        point_ind_on_edge2 = sort(unique(vec(mesh.edge[mesh.edge_marker.==pair[2],:]')))

        f = edge_trafo[i]

        point_edge1 = mesh.point[point_ind_on_edge1,:]
        # Apply f to edge1, gives permuted edge2
        point_edge2_on_edge1 = f(mesh.point[point_ind_on_edge2,:])

        # Find the indices of mapped points on edge2 in edge1, then change the
        # index of points on edge2 to get identified with edge1
        for j=1:size(point_edge2_on_edge1,1)
            ind_p2_in_p1 = closest_index(point_edge1, point_edge2_on_edge1[j,:])
            
            ind_of_point_found_in_global_vec = find(map(x->x[1]==point_ind_on_edge2[ind_p2_in_p1],index_map))[]

            # index_map[ind_of_point_found_in_global_vec][end+1] = point_ind_on_edge1[j]
            push!(index_map[ind_of_point_found_in_global_vec],point_ind_on_edge1[j])
        end
    end
    index_map = map(x->unique(sort(x)), index_map)

    for k=1:length(index_map)
        while length(index_map[k])>1
            ind_2b_replaced = index_map[k][end]
            ind_new = index_map[k][1]
            
            map(x-> (x[x.==ind_2b_replaced]=ind_new), index_map);
           
            # Keep only unique values
            index_map = map(x->unique(x), index_map)
        end        
    end

    index_map = vcat(index_map...)

    index_map_all = collect(1:mesh.n_point)
    index_map_all[ind_point_boundary] = index_map

    # Squeze dof indices
    d = setdiff(1:maximum(index_map_all), index_map_all)
    while length(d)>0
        index_map_all[index_map_all.>d[1]] -= 1
        d = setdiff(1:maximum(index_map_all), index_map_all)
    end

    return index_map_all
end

# find closest point in array
function closest_index(P :: Array{Float64,2}, p :: Array{Float64,1})

    ibest = start(eachindex(P[:,1]))
    dxbest = sum(abs.(P[ibest,:] - p))

    for ind in eachindex(P[:,1])
        dx = sum(abs.(P[ind,:] - p))
        if dx < dxbest
            dxbest = dx
            ibest = ind
        end
    end
    
    return ibest
end
# ----------------------------------------