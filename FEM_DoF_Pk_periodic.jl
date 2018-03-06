struct Dof_Pk_periodic{FEM_order} <: AbstractDof
    
    # ----------------------------------------
    # General infos
    n_order :: Int64
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
    n_edge_dirichlet :: Int64
    n_edge_neumann :: Int64
    n_edge_per_elem :: Int64

    ind_edge_boundary :: Array{Int64,1}
    ind_edge_interior :: Array{Int64,1}
    ind_edge_dirichlet :: Array{Int64,1}
    ind_edge_neumann :: Array{Int64,1}
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
    
    T_ref2cell :: Array{Array{Float64,2},1}
    T_cell2ref :: Array{Array{Float64,2},1}
    
end # end type



function Dof_Pk_periodic(mesh :: Mesh.TriangleMesh.TriMesh,
                            problem :: Problem.AbstractProblem,
                            periodicityInfo :: Mesh.PeriodicityInfo,
                            n_order :: Int)
    
    FEM_info = "DoF object for ---   periodic  --- Pk-Lagrange FEM of order $(n_order)."
    # ----------------------------------------

    if n_order==1
        # -------------------------------------------------------------------------------------------------

        # ----------------------------------------
        # Translation between topology and geometry
        index_map_dof2mesh = identify_points(mesh, periodicityInfo)
        index_map_mesh2dof = indexin(unique(index_map_dof2mesh), index_map_dof2mesh)
        # ----------------------------------------


        # ----------------------------------------
        # Node infos
        n_node = mesh.n_point

        label_edges_with_no_partner = setdiff(sort(unique(mesh.edge_marker[mesh.edge_marker.!=0])), 
                                                    sort(unique(periodicityInfo.edge_marker_pair[:])))
        label_edges_with_partner = sort(unique(periodicityInfo.edge_marker_pair[:]))

        ind_points_with_no_partner = unique(cat(1, [unique(vec(mesh.edge[:,mesh.edge_marker.==label])) for label in label_edges_with_no_partner]...))
        ind_points_with_partner = unique(cat(1, [unique(vec(mesh.edge[:,mesh.edge_marker.==label]))  for label in label_edges_with_partner]...))

        ind_node_boundary = ind_points_with_no_partner
        ind_node_interior = setdiff(1:n_node, ind_points_with_no_partner)

        ind_node_dirichlet = sort(unique(cat(1, [unique(mesh.edge[:,mesh.edge_marker.==marker]) for marker in problem.marker_dirichlet_edge]...)))
        ind_node_neumann = sort(unique(cat(1, [unique(mesh.edge[:,mesh.edge_marker.==marker]) for marker in problem.marker_neumann_edge]...)))
        ind_node_non_dirichlet = setdiff(1:n_node, ind_node_dirichlet)

        n_node_dirichlet = length(ind_node_dirichlet)
        n_node_neumann = length(ind_node_neumann)
        n_node_boundary = n_node_dirichlet + n_node_neumann
        n_node_interior = n_node - n_node_boundary
        n_node_per_edge = 2
        n_node_per_elem = 3
        # ----------------------------------------
        
        
        # ----------------------------------------
        # Edge infos
        n_edge = mesh.n_edge

        ind_edge_boundary  = unique(cat(1, [find(mesh.edge_marker.==label) for label in label_edges_with_no_partner]...))
        ind_edge_dirichlet  = unique(cat(1, [find(mesh.edge_marker.==label) for label in problem.marker_dirichlet_edge]...))
        ind_edge_neumann  = unique(cat(1, [find(mesh.edge_marker.==label) for label in problem.marker_neumann_edge]...))
        ind_edge_interior = setdiff(1:mesh.n_edge, ind_edge_boundary)

        n_edge_boundary = length(ind_edge_boundary)
        n_edge_interior = length(ind_edge_interior)
        n_edge_dirichlet = length(ind_edge_dirichlet)
        n_edge_neumann = length(ind_edge_neumann)
        n_edge_per_elem = 3
        # ----------------------------------------
        
        
        # ----------------------------------------
        # Element infos
        n_elem = mesh.n_cell
        
        # All elements with at least two boundary nodes
        ind_elem_boundary = find([length(setdiff(mesh.cell[:,i], ind_node_boundary)) for i=1:n_elem].==1)
        ind_elem_interior = setdiff(1:n_elem, ind_elem_boundary)

        n_elem_boundary = length(ind_elem_boundary)
        n_elem_interior = length(ind_elem_interior)
        # ----------------------------------------
        
        
        # ----------------------------------------
        # Topology info
        is_periodic = true

        n_true_dof = length(unique(index_map_dof2mesh))
        # ----------------------------------------

        # ----------------------------------------
        # Used for building the system matrices
        ind_cell = index_map_dof2mesh[Mesh.get_cell(mesh, 1:mesh.n_cell)]
        ind = vec(ind_cell[[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3],:])
        ind_test = vec(ind_cell[[1;2;3;1;2;3;1;2;3],:])
        ind_lin = sub2ind((n_true_dof,n_true_dof), ind_test, ind)
        
        T_ref2cell = [zeros(2, 3) for i=1:mesh.n_cell]
        T_cell2ref = [zeros(2, 3) for i=1:mesh.n_cell]
        for i=1:n_elem
            # Extended transformation matrices for mapping from an to the
            # reference cell K = [(0,0), (1,0), (0,1)]
            P = Mesh.get_point(mesh, Mesh.get_cell(mesh,[i]))

            T_ref2cell[i] = [P[:,2]-P[:,1]  P[:,3]-P[:,1] P[:,1]]
            T_cell2ref[i] = (eye(2)/T_ref2cell[i][:,1:2]) * [   eye(2)  mesh.point[:,mesh.cell[1,i]]   ]
        end
        # ----------------------------------------
        # -------------------------------------------------------------------------------------------------
            
    elseif n_order==2
        # -------------------------------------------------------------------------------------------------
        error("Order 2 DoF not implemented yet (periodic setting).")
        # -------------------------------------------------------------------------------------------------

    elseif n_order==3
        # -------------------------------------------------------------------------------------------------
        error("Order 3 DoF not implemented yet (periodic setting).")
        # -------------------------------------------------------------------------------------------------

    end # end if order
    
    return Dof_Pk_periodic{n_order}(n_order,
                                        FEM_info,
                                        n_node,
                                        n_node_boundary,
                                        n_node_interior,
                                        n_node_dirichlet,
                                        n_node_neumann,
                                        n_node_per_edge,
                                        n_node_per_elem,
                                        ind_node_boundary,
                                        ind_node_interior,
                                        ind_node_dirichlet,
                                        ind_node_non_dirichlet,
                                        ind_node_neumann,
                                        n_edge,
                                        n_edge_boundary,
                                        n_edge_interior,
                                        n_edge_dirichlet,
                                        n_edge_neumann,
                                        n_edge_per_elem,
                                        ind_edge_boundary,
                                        ind_edge_interior,
                                        ind_edge_dirichlet,
                                        ind_edge_neumann,
                                        n_elem,
                                        n_elem_boundary,
                                        n_elem_interior,
                                        ind_elem_boundary,
                                        ind_elem_interior,
                                        is_periodic,
                                        n_true_dof,
                                        index_map_dof2mesh,
                                        index_map_mesh2dof,
                                        ind,
                                        ind_test,
                                        ind_lin,
                                        T_ref2cell,
                                        T_cell2ref)
end # end constructor



# -------------------------------------------------------------------------------------------------------
# -----------------------------   Functions on Dof_Pk_periodic   -----------------------------
# -------------------------------------------------------------------------------------------------------

# ----------------------------------------
"""
    map_ind_dof2mesh(dof :: Dof_Pk_periodic{1}, ind_dof :: Array{Int})

    Map indices in the actual mesh to indices in terms of degrees of freedom.
    
"""
function map_ind_dof2mesh(dof :: Dof_Pk_periodic{1}, ind_dof :: Array{Int})

    ind_mesh = [find(index_map_dof2mesh.==ind) for ind in ind_dof]

    return ind_dof
end


"""
    map_vec_dof2mesh(dof :: Dof_Pk_periodic{1}, vec_dof :: Array{Float64,1})

    Map a vector in terms of the actual mesh to a vector in terms of degrees
    of freedom.
    
    
"""
function map_vec_dof2mesh(dof :: Dof_Pk_periodic{1}, vec_dof :: Array{Float64,1})

    # expand a dof-vector into a mesh-vector (only periodic boundaries)
    vec_mesh = vec_dof[dof.index_map_dof2mesh]

    return vec_mesh
end


"""
    map_vec_dof2mesh(dof :: Dof_Pk_periodic{1}, vec_dof :: Array{Float64,2})

    Map a vector in terms of the actual mesh to a vector in terms of degrees
    of freedom.
    
    
"""
function map_vec_dof2mesh(dof :: Dof_Pk_periodic{1}, vec_dof :: Array{Float64,2})

    # expand a dof-vector into a mesh-vector (only periodic boundaries)
    vec_mesh = vec_dof[dof.index_map_dof2mesh,:]

    return vec_mesh
end
# ----------------------------------------


# ----------------------------------------
"""
    map_ind_mesh2dof(dof :: Dof_Pk_periodic{1}, ind_mesh :: Array{Int})

    Map indices in terms of degrees of freedom to indices on the actual
    mesh.

"""
function map_ind_mesh2dof(dof :: Dof_Pk_periodic{1}, ind_mesh :: Array{Int})

    ind_dof = dof.index_map_dof2mesh[ind_mesh]

    return ind_dof
end


"""
    map_vec_mesh2dof(dof :: Dof_Pk_periodic{1}, vec_mesh :: Array{Float64,1})

    Map a vector in terms of degrees of freedom to a vector on the actual
    mesh.

"""
function map_vec_mesh2dof(dof :: Dof_Pk_periodic{1}, vec_mesh :: Array{Float64,1})

    # reduce a mesh-vector to a dof-vector (only periodic boundaries)
    vec_dof = vec_mesh[dof.index_map_mesh2dof]

    return vec_dof
end


"""
    map_vec_mesh2dof(dof :: Dof_Pk_periodic{1}, vec_mesh :: Array{Float64,2})

    Map a vector in terms of degrees of freedom to a vector on the actual
    mesh.

"""
function map_vec_mesh2dof(dof :: Dof_Pk_periodic{1}, vec_mesh :: Array{Float64,2})

    # reduce a mesh-vector to a dof-vector (only periodic boundaries)
    vec_dof = vec_mesh[dof.index_map_mesh2dof,:]

    return vec_dof
end
# ----------------------------------------


# ----------------------------------------
function get_dof_elem(dof :: Dof_Pk_periodic{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_c :: Array{Int64,1})
    
    # Put a filter before mesh indices to translate to dof
    # indices
    
    return map_ind_mesh2dof(dof, Mesh.get_cell(mesh, ind_c))
end # end function

function get_dof_elem(dof :: Dof_Pk_periodic{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_c :: UnitRange{Int64})
    
    # Put a filter before mesh indices to translate to dof
    # indices
    
    return get_dof_elem(dof, mesh, collect(ind_c))
end # end function
# ----------------------------------------


# ----------------------------------------
function get_dof_edge(dof :: Dof_Pk_periodic{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_e :: Array{Int64,1})    
    
    # Put a filter before mesh indices to translate to dof
    # indices
    
    return map_ind_mesh2dof(dof, Mesh.get_edge(mesh, ind_e))
end # end function

function get_dof_edge(dof :: Dof_Pk_periodic{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_e :: UnitRange{Int64})
    
    # Put a filter before mesh indices to translate to dof
    # indices
    
    return get_dof_edge(dof, mesh, collect(ind_e))
end # end function
# ----------------------------------------


# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------


# ----------------------------------------
# ----------------------------------------
function identify_points(mesh :: Mesh.TriangleMesh.TriMesh, periodicityInfo :: Mesh.PeriodicityInfo)

    edge_marker_pair = periodicityInfo.edge_marker_pair
    edge_trafo =  periodicityInfo.f_edge_to_edge

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

        point_ind_on_edge1 = sort(unique(vec(mesh.edge[:,mesh.edge_marker.==pair[1]])))
        point_ind_on_edge2 = sort(unique(vec(mesh.edge[:,mesh.edge_marker.==pair[2]])))

        f = edge_trafo[i]

        point_edge1 = mesh.point[:,point_ind_on_edge1]
        # Apply f to edge2, gives permuted edge1
        point_edge2_on_edge1 = f(mesh.point[:,point_ind_on_edge2])

        # Find the indices of mapped points of edge2 in edge1, then change the
        # index of points on edge2 to get identified with edge1
        for j=1:size(point_edge2_on_edge1,2)
            ind_p2_in_p1 = closest_index(point_edge1, point_edge2_on_edge1[:,j])
            # This means point_ind_on_edge2[j] -> point_ind_on_edge1[ind_p2_in_p1]

            # Index of point_ind_on_edge2[j] must now be found in index_map
            ind_of_point_found_in_global_vec = find(map(x->x[1]==point_ind_on_edge2[j],index_map))[]
            
            # Then we change add to
            # index_map[ind_of_point_found_in_global_vec] the identified
            # corresponding point
            push!(index_map[ind_of_point_found_in_global_vec],point_ind_on_edge1[ind_p2_in_p1])
        end
    end
    
    # sort the list of points that are identified with each other
    index_map = map(x->unique(sort(x)), index_map)

    # For each point indentified with a point of a lower index chnage the
    # identifier of that point to the lower index
    for k=1:length(index_map)
        while length(index_map[k])>1
            ind_2b_replaced = index_map[k][end]
            ind_new = index_map[k][1]
            
            map(x-> (x[x.==ind_2b_replaced]=ind_new), index_map);
           
            # Keep only unique values
            index_map = map(x->unique(x), index_map)
        end        
    end

    # convert this to a column vector
    index_map = vcat(index_map...)

    # The indexing has gaps. We need to reduce indices as long as thera are
    # gaps.
    index_map_all = collect(1:mesh.n_point)
    index_map_all[ind_point_boundary] = index_map

    # This is the actual reduction.
    d = setdiff(1:maximum(index_map_all), index_map_all)
    while length(d)>0
        index_map_all[index_map_all.>d[1]] -= 1
        d = setdiff(1:maximum(index_map_all), index_map_all)
    end

    return index_map_all
end

# find closest point in array
function closest_index(P :: Array{Float64,2}, p :: Array{Float64,1})

    ibest = start(eachindex(P[1,:]))
    dxbest = sum(abs.(P[:,ibest] - p))

    for ind in eachindex(P[1,:])
        dx = sum(abs.(P[:,ind] - p))
        if dx < dxbest
            dxbest = dx
            ibest = ind
        end
    end
    
    return ibest
end
# ----------------------------------------
# ----------------------------------------