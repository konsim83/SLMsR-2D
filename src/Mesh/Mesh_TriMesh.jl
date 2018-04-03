# ---------------------------------------------------------------------------------------------
# -------------------   Additional constructor for TriangleMesh.TriMesh   -------------------
# ---------------------------------------------------------------------------------------------
function TriMesh(mesh :: TriangleMesh.TriMesh, point :: Array{Float64,2}, mesh_info :: String)

    voronoi = TriangleMesh.VoronoiDiagram("No Voronoi diagram", 
                                            0, Array{Float64,2}(2,0),
                                            0, Array{Float64,2}(0,0),
                                            0, Array{Int,2}(2,0),
                                            Array{Int,2}(2,0))
    n_hole = 0
    hole = Array{Float64,2}(2,0)

    return TriangleMesh.TriMesh(mesh_info,
                                mesh.n_point, point,
                                mesh.n_point_marker, mesh.point_marker,
                                mesh.n_point_attribute, mesh.point_attribute,
                                mesh.n_cell, mesh.cell, mesh.cell_neighbor,
                                mesh.n_edge, mesh.edge, 
                                mesh.n_edge_marker, mesh.edge_marker,
                                mesh.n_segment, mesh.segment,
                                mesh.n_segment_marker, mesh.segment_marker,
                                n_hole, hole,
                                voronoi)
end # end constructor
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------
# -----------------------------   Functions on TriMesh   -----------------------------
# --------------------------------------------------------------------------------------------

# -------   Two versions of get_cell   -------
function get_cell(m :: TriangleMesh.TriMesh, ind_c :: Array{Int64,1} )
    
    return m.cell[:,ind_c]
end # end get_cell


function get_cell(m :: TriangleMesh.TriMesh, ind_c :: UnitRange{Int64} )
    
    return get_cell(m :: TriangleMesh.TriMesh, collect(ind_c) )
end # end get_cell
# -------   Two versions of get_cell   -------



# -------   Two versions of get_edge   -------
function get_edge( m :: TriangleMesh.TriMesh, ind_e :: Array{Int64,1} )
    
    return m.edge[:,ind_e]
end # end get_edge


function get_edge(m :: TriangleMesh.TriMesh, ind_e :: UnitRange{Int64} )
    
    return get_edge(m, collect(ind_e))
end # end get_edge
# -------   Two versions of get_edge   -------



# -------   Three versions of get_point   -------
function get_point(m :: TriangleMesh.TriMesh, ind_p :: Array{Int64,1})
    """ Version for type Vector """
            
    return m.point[:,ind_p]
end # end get_point


function get_point(m :: TriangleMesh.TriMesh, ind_p :: UnitRange{Int64} )
    """ Version for type UnitRange """
    
    return get_point(m, collect(ind_p))
end # end get_point


function get_point(m :: TriangleMesh.TriMesh, ind_p :: Array{Int64,2})
    """

    Version for connectivity list input. The dimensions are
    (n_point_per_cell, 2, n_cell)
    
    """
    P = Array{Float64,3}(2, 3, size(ind_p,2))

    P[:,:,:] = m.point[:,ind_p]

    return P
end # end get_point
# -------   Three versions of get_point   -------

# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
