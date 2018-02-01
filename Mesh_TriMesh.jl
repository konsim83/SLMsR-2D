# Actual Julia mesh type
mutable struct TriMesh
    point :: Array{Float64, 2}
    point_marker :: Array{Int64, 2}
    point_attribute :: Array{Float64, 2}
    n_point :: Int64

    cell :: Array{Int64, 2}
    cell_neighbor :: Array{Int64, 2}
    n_cell :: Int64

    edge :: Array{Int64, 2}
    edge_marker :: Array{Int64, 1}
    n_edge :: Int64

    segment :: Array{Int64, 2}
    segment_marker :: Array{Int64, 1}
    n_segment :: Int64

    T_ref2cell :: Array{Float64,3}
    T_cell2ref :: Array{Float64,3}
    
    mesh_info :: String

    # ---------------------------------------------------------------------------------------------
    function TriMesh(mesh :: TriMesh, point :: Array{Float64,2}, info :: String)
        this = new()

        this.mesh_info = info
        
        # Points
        this.n_point = size(point,1)
        this.point = point
        this.point_marker = mesh.point_marker
        this.point_attribute = mesh.point_attribute
                
        #Triangles
        this.n_cell = mesh.n_cell
        this.cell = mesh.cell

        this.T_ref2cell = zeros(2, 3, mesh.n_cell)
        this.T_cell2ref = zeros(2, 3, mesh.n_cell)
        this.cell_neighbor = mesh.cell_neighbor

        for i=1:this.n_cell            
            # Extended transformation matrices for mapping from an to the reference cell K = [(0,0), (1,0), (0,1)]
            this.T_ref2cell[:,:,i] = [   this.point[this.cell[i,2],:]-this.point[this.cell[i,1],:]   this.point[this.cell[i,3],:]-this.point[this.cell[i,1],:]   this.point[this.cell[i,1],:]   ]
            this.T_cell2ref[:,:,i] = (eye(2)/this.T_ref2cell[:,1:2,i]) * [   eye(2)  this.point[this.cell[i,1],:]   ]
        end
        
        # Edges
        this.n_edge = mesh.n_edge
        this.edge = mesh.edge
        this.edge_marker = mesh.edge_marker

        # Segments
        this.n_segment = mesh.n_segment
        this.segment = mesh.segment
        this.segment_marker = mesh.segment_marker

        return this
    end # end constructor
    # ---------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------
    function TriMesh(mesh :: TriangleMesh.TriMesh)
        this = new()

        this.mesh_info = mesh.mesh_info
        
        # Points
        this.n_point = size(mesh.point,1)
        this.point = mesh.point
        this.point_marker = mesh.point_marker
        this.point_attribute = mesh.point_attribute
                
        #Triangles
        this.n_cell = mesh.n_cell
        this.cell = mesh.cell

        this.T_ref2cell = zeros(2, 3, mesh.n_cell)
        this.T_cell2ref = zeros(2, 3, mesh.n_cell)
        this.cell_neighbor = mesh.cell_neighbor

        for i=1:this.n_cell            
            # Extended transformation matrices for mapping from an to the reference cell K = [(0,0), (1,0), (0,1)]
            this.T_ref2cell[:,:,i] = [   this.point[this.cell[i,2],:]-this.point[this.cell[i,1],:]   this.point[this.cell[i,3],:]-this.point[this.cell[i,1],:]   this.point[this.cell[i,1],:]   ]
            this.T_cell2ref[:,:,i] = (eye(2)/this.T_ref2cell[:,1:2,i]) * [   eye(2)  this.point[this.cell[i,1],:]   ]
        end
        
        # Edges
        this.n_edge = mesh.n_edge
        this.edge = mesh.edge
        this.edge_marker = mesh.edge_marker

        # Segments
        this.n_segment = mesh.n_segment
        this.segment = mesh.segment
        this.segment_marker = mesh.segment_marker

        return this
    end # end constructor
    # ---------------------------------------------------------------------------------------------
end # end type



# --------------------------------------------------------------------------------------------
# -----------------------------   Functions on TriMesh   -----------------------------
# --------------------------------------------------------------------------------------------

# -------   Two versions of get_cell   -------
function get_cell( m :: TriMesh, ind_c :: UnitRange{Int64} )
    
    return get_cell( m :: TriMesh, collect(ind_c) )
end # end get_cell


function get_cell( m :: TriMesh, ind_c :: Array{Int64,1} )
    
    return m.cell[ind_c,:]
end # end get_cell
# -------   Two versions of get_cell   -------



# -------   Two versions of get_edge   -------
function get_edge(m :: TriMesh, ind_e :: UnitRange{Int64} )
    
    return get_edge(m, collect(ind_e))
end # end get_edge


function get_edge( m :: TriMesh, ind_e :: Array{Int64,1} )
    
    return m.edge[ind_e,:]
end # end get_edge
# -------   Two versions of get_edge   -------



# -------   Three versions of get_point   -------
function get_point( m :: TriMesh, ind_p :: UnitRange{Int64} )
    """ Version for type UnitRange """
    
    return get_point(m, collect(ind_p))
end # end get_point


function get_point( m :: TriMesh, ind_p :: Array{Int64,1})
    """ Version for type Vector """
            
    return m.point[ind_p,:]
end # end get_point


function get_point( m :: TriMesh, ind_p :: Array{Int64,2})
    """

    Version for connectivity list input. The dimensions are
    (n_point_per_cell, 2, n_cell)
    
    """
    P = Array{Float64,3}(3, 2, size(ind_p,1))

    P[:,1,:] = m.point[ind_p',1]
    P[:,2,:] = m.point[ind_p',2]
    
    return P
end # end get_point
# -------   Three versions of get_point   -------

# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------
