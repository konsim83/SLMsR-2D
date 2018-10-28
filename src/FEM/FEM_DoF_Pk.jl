struct Dof_Pk{FEM_order} <: AbstractDof
    
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
    # Topology info
    is_periodic :: Bool
    n_true_dof :: Int64

    map_vec_ind_mesh2dof :: Array{Int64,1}
    # ----------------------------------------

    ind :: Array{Int64,1}
    ind_test :: Array{Int64,1}
    ind_lin :: Array{Int64,1}

    T_ref2cell :: Array{Array{Float64,2},1}
    T_cell2ref :: Array{Array{Float64,2},1}

end # end type



function Dof_Pk(mesh :: Mesh.TriangleMesh.TriMesh,
                problem :: Problem.AbstractProblem,
                n_order :: Int)
    
    FEM_info = "DoF object for ---   non-periodic   --- Pk-Lagrange FEM of order $(n_order)."
    

    if n_order==1
        # -------------------------------------------------------------------------------------------------
        
        # ----------------------------------------
        # Node infos
        n_node = mesh.n_point

        ind_node_boundary = findall(mesh.point_marker.!=0)
        ind_node_interior = findall(mesh.point_marker.==0)
        ind_node_dirichlet  = sort(unique(cat(dims=1, [unique(mesh.edge[:,mesh.edge_marker.==marker]) for marker in problem.marker_dirichlet_edge]...)))
        ind_node_non_dirichlet = setdiff(1:n_node, ind_node_dirichlet)
        ind_node_neumann  = setdiff(ind_node_non_dirichlet, ind_node_interior)

        n_node_boundary = sum(mesh.point_marker.!=0)
        n_node_interior = sum(mesh.point_marker.==0)
        n_node_dirichlet = length(ind_node_dirichlet)
        n_node_neumann = n_node_boundary - n_node_dirichlet
        n_node_per_edge = 2
        n_node_per_elem = 3
        # ----------------------------------------
        
        
        # ----------------------------------------
        # Edge infos
        n_edge = mesh.n_edge

        ind_edge_boundary = findall(mesh.edge_marker.!=0)
        ind_edge_interior = findall(mesh.edge_marker.==0)
        ind_edge_dirichlet = sort(cat(dims=1, [findall(mesh.edge_marker.==marker) for marker in problem.marker_dirichlet_edge]...))
        ind_edge_neumann = sort(cat(dims=1, [findall(mesh.edge_marker.==marker) for marker in problem.marker_neumann_edge]...))

        
        n_edge_boundary = length(ind_edge_boundary)
        n_edge_interior = length(ind_edge_interior)
        n_edge_dirichlet = length(ind_edge_dirichlet)
        n_edge_neumann = length(ind_edge_neumann)
        n_edge_per_elem = 3
        # ----------------------------------------
        
        
        # ----------------------------------------
        # Element infos
        n_elem = mesh.n_cell
        
        ind_elem_boundary = findall(sum(mesh.cell_neighbor.==0,dims=1).!=0)
        ind_elem_interior = findall(sum(mesh.cell_neighbor.==0,dims=1).==0)

        n_elem_boundary = sum(sum(mesh.cell_neighbor.==0,dims=1).!=0)
        n_elem_interior = sum(sum(mesh.cell_neighbor.==0,dims=1).==0)
        # ----------------------------------------

        
        # ----------------------------------------
        # Topology info
        is_periodic = false
        n_true_dof = mesh.n_point

        # No nodes need to be identified with each other
        map_vec_ind_mesh2dof = collect(1:n_node)
        # ----------------------------------------


        # ----------------------------------------
        # Used for building the system matrices
        ind_cell = Mesh.get_cell(mesh, 1:mesh.n_cell)
        ind = vec(ind_cell[[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3],:])
        ind_test = vec(ind_cell[[1;2;3;1;2;3;1;2;3],:])

        # Convert set of cartesian coordinates to linear indices
        ic = CartesianIndices((n_true_dof,n_true_dof))
        il = LinearIndices((n_true_dof,n_true_dof))
        indCart = [ic[ind_test[k],ind[k]] for k in 1:length(ind_test)]
        ind_lin = il[indCart]

        T_ref2cell = [zeros(2, 3) for i=1:mesh.n_cell]
        T_cell2ref = [zeros(2, 3) for i=1:mesh.n_cell]
        for i=1:n_elem
            # Extended transformation matrices for mapping from an to the
            # reference cell K = [(0,0), (1,0), (0,1)]
            P = Mesh.get_point(mesh, Mesh.get_cell(mesh,[i]))

            T_ref2cell[i] = [P[:,2]-P[:,1]  P[:,3]-P[:,1] P[:,1]]
            T_cell2ref[i] = (I/T_ref2cell[i][:,1:2]) * [   I  mesh.point[:,mesh.cell[1,i]]   ]
        end
        # ----------------------------------------
        # -------------------------------------------------------------------------------------------------

    elseif n_order==2
        # -------------------------------------------------------------------------------------------------
        error("Order 2 DoF not implemented yet.")
        # -------------------------------------------------------------------------------------------------
    
    elseif n_order==3
        # -------------------------------------------------------------------------------------------------
        error("Order 3 DoF not implemented yet.")
        # -------------------------------------------------------------------------------------------------

    end # end if order

    return Dof_Pk{n_order}(n_order,
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
                            map_vec_ind_mesh2dof,
                            ind,
                            ind_test,
                            ind_lin,
                            T_ref2cell,
                            T_cell2ref)
end



# -------------------------------------------------------------------------------------------
# -----------------------------   Functions on Dof_Pk   -----------------------------
# -------------------------------------------------------------------------------------------

# ----------------------------------------
"""
    map_ind_dof2mesh(dof :: Dof_Pk{1}, ind_dof :: Array{Int})

    Map indices in the actual mesh to indices in terms of degrees of freedom.
    
"""
function map_ind_dof2mesh(dof :: Dof_Pk{1}, ind_dof :: Array{Int})

    return ind_dof
end


"""
    map_vec_dof2mesh(dof :: Dof_Pk{1}, vec_dof :: Array{Float64,1})

    Map a vector in terms of the actual mesh to a vector in terms of degrees
    of freedom.
    
    
"""
function map_vec_dof2mesh(dof :: Dof_Pk{1}, vec_dof :: Array{Float64,1})

    return vec_dof
end


"""
    map_vec_dof2mesh(dof :: Dof_Pk{1}, vec_dof :: Array{Float64,2})

    Map a vector in terms of the actual mesh to a vector in terms of degrees
    of freedom.
    
    
"""
function map_vec_dof2mesh(dof :: Dof_Pk{1}, vec_dof :: Array{Float64,2})

    return vec_dof
end
# ----------------------------------------


# ----------------------------------------
"""
    map_ind_mesh2dof(dof :: Dof_Pk{1}, ind_mesh :: Array{Int})

    Map indices in terms of degrees of freedom to indices on the actual
    mesh.

"""
function map_ind_mesh2dof(dof :: Dof_Pk{1}, ind_mesh :: Array{Int})

    return ind_mesh
end


"""
    map_vec_mesh2dof(dof :: Dof_Pk{1}, vec_mesh :: Array{Float64,1})

    Map a vector in terms of degrees of freedom to a vector on the actual
    mesh.

"""
function map_vec_mesh2dof(dof :: Dof_Pk{1}, vec_mesh :: Array{Float64,1})

    return vec_mesh
end


"""
    map_vec_mesh2dof(dof :: Dof_Pk{1}, vec_mesh :: Array{Float64,2})

    Map a vector in terms of degrees of freedom to a vector on the actual
    mesh.

"""
function map_vec_mesh2dof(dof :: Dof_Pk{1}, vec_mesh :: Array{Float64,2})

    return vec_mesh
end
# ----------------------------------------


# ----------------------------------------
function get_dof_elem(dof :: Dof_Pk{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_c :: Array{Int64,1})
    
    return Mesh.get_cell(mesh, ind_c)
end # end function

function get_dof_elem(dof :: Dof_Pk{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_c :: UnitRange{Int64})
    
    return Mesh.get_cell(mesh, ind_c)
end # end function
# ----------------------------------------


# ----------------------------------------
function get_dof_edge(dof :: Dof_Pk{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_e :: Array{Int64,1})
    
    return Mesh.get_edge(mesh, ind_e)
end # end function

function get_dof_edge(dof :: Dof_Pk{1}, mesh :: Mesh.TriangleMesh.TriMesh, ind_e :: UnitRange{Int64})
    
    return Mesh.get_edge(mesh, ind_e)
end # end function
# ----------------------------------------

# -------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------