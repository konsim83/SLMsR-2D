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

    get_point :: Function
    get_edge :: Function
    get_cell :: Function

    mesh_info :: String
    
    # ---------------------------------------------------------------------------------------------
    function TriMesh(mesh :: Triangle_mesh_C, info :: String)
        this = new()

        this.mesh_info = info
        
        # Points
        this.n_point = mesh.numberofpoints
        this.point = zeros(mesh.numberofpoints, 2)
        this.point_marker = zeros(mesh.numberofpoints)
        for i=1:this.n_point
            this.point[i,1] = unsafe_load(mesh.pointlist, 2*i-1);
            this.point[i,2] = unsafe_load(mesh.pointlist, 2*i);
            this.point_marker[i] = unsafe_load(mesh.pointmarkerlist, i);
        end
                
        #Triangles
        this.n_cell = mesh.numberoftriangles
        this.cell = zeros(mesh.numberoftriangles, 3)
        this.T_ref2cell = zeros(2, 3, mesh.numberoftriangles)
        this.T_cell2ref = zeros(2, 3, mesh.numberoftriangles)
        this.cell_neighbor = zeros(mesh.numberoftriangles, 3)
        for i=1:this.n_cell
            this.cell[i,1] = unsafe_load(mesh.trianglelist, 3*i-2);
            this.cell[i,2] = unsafe_load(mesh.trianglelist, 3*i-1);
            this.cell[i,3] = unsafe_load(mesh.trianglelist, 3*i);
            
            this.cell_neighbor[i,1] = unsafe_load(mesh.neighborlist, 3*i-2);
            this.cell_neighbor[i,2] = unsafe_load(mesh.neighborlist, 3*i-1);
            this.cell_neighbor[i,3] = unsafe_load(mesh.neighborlist, 3*i);
            
            # Extended transformation matrices for mapping from an to the reference cell K = [(0,0), (1,0), (0,1)]
            this.T_ref2cell[:,:,i] = [   this.point[this.cell[i,2]+1,:]-this.point[this.cell[i,1]+1,:]   this.point[this.cell[i,3]+1,:]-this.point[this.cell[i,1]+1,:]   this.point[this.cell[i,1]+1,:]   ]
            this.T_cell2ref[:,:,i] = (eye(2)/this.T_ref2cell[:,1:2,i]) * [   eye(2)  this.point[this.cell[i,1]+1,:]   ]
        end
        this.cell_neighbor = this.cell_neighbor + 1;
        this.cell = this.cell + 1;
        
        # Edges
        this.n_edge = mesh.numberofedges
        this.edge = zeros(mesh.numberofedges, 2)
        this.edge_marker = zeros(mesh.numberofedges)
        for i=1:this.n_edge
            this.edge[i,1] = unsafe_load(mesh.edgelist, 2*i-1);
            this.edge[i,2] = unsafe_load(mesh.edgelist, 2*i);
            this.edge_marker[i] = unsafe_load(mesh.edgemarkerlist, i);
        end
        this.edge = this.edge + 1;

        # Segments
        this.n_segment = mesh.numberofsegments
        this.segment = zeros(mesh.numberofsegments, 2)
        this.segment_marker = zeros(mesh.numberofedges)
        for i=1:this.n_segment
            this.segment[i,1] = unsafe_load(mesh.segmentlist, 2*i-1);
            this.segment[i,2] = unsafe_load(mesh.segmentlist, 2*i);
            this.segment_marker[i] = unsafe_load(mesh.segmentmarkerlist, i);
        end
        this.segment = this.segment + 1;
        

        this.get_cell = function( ind_c )
            ind_c = vec(collect(ind_c))
            
            return this.cell[ind_c,:]
        end # end get_cell

        this.get_edge = function( ind_e )
            ind_e = vec(collect(ind_e))
            
            return this.edge[ind_e,:]
        end # end get_edge
            
        this.get_point = function( ind_p )
            
            return this.point[ind_p,:]
        end # end get_point
        
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

        this.get_cell = function( ind_c )
            ind_c = vec(collect(ind_c))
            
            return this.cell[ind_c,:]
        end # end get_cell

        this.get_edge = function( ind_e )
            ind_e = vec(collect(ind_e))
            
            return this.edge[ind_e,:]
        end # end get_edge
            
        this.get_point = function( ind_p )
            
            return this.point[ind_p,:]
        end # end get_point
        
        return this
    end # end constructor
    # ---------------------------------------------------------------------------------------------
end # end type
