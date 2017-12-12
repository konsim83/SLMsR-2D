# Actual Julia mesh type
type TriMesh
    point :: Array{Float64, 2}
    point_marker :: Array{Int64, 1}
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
    function TriMesh(mesh :: Triangle_mesh_C, info :: String)
        this = new()

        this.mesh_info = string("Triangular mesh of " , info)


        # ------------------ This function has two bottlenecks:
        # 1. Menory layout of C vs. Fortran while loading the pointers
        # 2. The loop is slow
        # ------------------
        

        # Points
        this.n_point = mesh.numberofpoints
        
        this.point = Array{Float64,2}(this.n_point, 2)
        this.point = unsafe_wrap(Array, mesh.pointlist, (2,this.n_point))'

        this.point_marker = Array{Int64,1}(this.n_point)
        this.point_marker = unsafe_wrap(Array, mesh.pointmarkerlist, this.n_point)

        
        #Triangles
        this.n_cell = mesh.numberoftriangles
        this.T_ref2cell = zeros(2, 3, mesh.numberoftriangles)
        this.T_cell2ref = zeros(2, 3, mesh.numberoftriangles)

        this.cell = Array{Int64,2}(this.n_cell, 3)
        this.cell = unsafe_wrap(Array, mesh.trianglelist, (3, this.n_cell))'

        this.cell_neighbor = Array{Int64,2}(this.n_cell, 3)
        this.cell_neighbor = unsafe_wrap(Array, mesh.neighborlist, (3,this.n_cell))'

        for i=1:this.n_cell
            # Extended transformation matrices for mapping from and to the reference cell K = [(0,0), (1,0), (0,1)]
            this.T_ref2cell[:,:,i] = [   this.point[this.cell[i,2]+1,:]-this.point[this.cell[i,1]+1,:]   this.point[this.cell[i,3]+1,:]-this.point[this.cell[i,1]+1,:]   this.point[this.cell[i,1]+1,:]   ]
            this.T_cell2ref[:,:,i] = (eye(2)/this.T_ref2cell[:,1:2,i]) * [   eye(2)  -this.point[this.cell[i,1]+1,:]   ]
        end
        
        this.cell_neighbor = this.cell_neighbor + 1;
        this.cell = this.cell + 1;

        
        # Edges
        this.n_edge = mesh.numberofedges

        this.edge = Array{Int64,2}(this.n_edge, 2)
        this.edge = unsafe_wrap(Array, mesh.edgelist, (2, this.n_edge))'
        
        this.edge_marker = Array{Int64,1}(this.n_edge)
        this.edge_marker = unsafe_wrap(Array, mesh.edgemarkerlist, this.n_edge)

        this.edge = this.edge + 1;


        # Segments
        this.n_segment = mesh.numberofsegments

        this.segment = Array{Int64,2}(this.n_segment, 2)
        this.segment = unsafe_wrap(Array, mesh.segmentlist, (2,this.n_segment))'

        this.segment_marker = Array{Int64,1}(this.n_segment)
        this.segment_marker = unsafe_wrap(Array, mesh.segmentmarkerlist, this.n_segment)

        this.segment = this.segment + 1;
        
        return this
    end # end constructor
    # ---------------------------------------------------------------------------------------------
            
    
    # ---------------------------------------------------------------------------------------------
    function TriMesh(mesh :: TriMesh, point :: Array{Float64,2}, info :: String)
        this = new()

        this.mesh_info = info
        
        # Points
        this.n_point = size(point,1)
        this.point = point
        this.point_marker = mesh.point_marker
                
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
