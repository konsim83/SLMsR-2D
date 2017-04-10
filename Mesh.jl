module Mesh

using PyPlot, WriteVTK

# ---------------------------------------------------------------------------------------------

# Type for reading C output
type Triangle_mesh_C
    pointlist :: Ptr{Cdouble}
    pointattributelist :: Ptr{Cdouble}
    pointmarkerlist :: Ptr{Cint}
    numberofpoints :: Cint
    numberofpointattributes :: Cint

    trianglelist :: Ptr{Cint}
    triangleattributelist :: Ptr{Cdouble}
    trianglearealist :: Ptr{Cdouble}
    neighborlist :: Ptr{Cint}
    numberoftriangles :: Cint
    numberofcorners :: Cint
    numberoftriangleattributes :: Cint

    segmentlist :: Ptr{Cint}
    segmentmarkerlist :: Ptr{Cint}
    numberofsegments :: Cint

    holelist :: Ptr{Cdouble}
    numberofholes :: Cint

    regionlist :: Ptr{Cdouble}
    numberofregions :: Cint

    edgelist :: Ptr{Cint}
    edgemarkerlist :: Ptr{Cint}
    normlist :: Ptr{Cdouble}
    numberofedges :: Cint
end # end type

# ---------------------------------------------------------------------------------------------

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
    n_edges :: Int64

    T_ref2cell :: Array{Float64,3}
    T_cell2ref :: Array{Float64,3}

    get_point :: Function
    get_cell :: Function
    
    function TriMesh(mesh :: Triangle_mesh_C)
        this = new()

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
        this.n_edges = mesh.numberofedges
        this.edge = zeros(mesh.numberofedges, 2)
        this.edge_marker = zeros(mesh.numberofedges)
        for i=1:this.n_edges
            this.edge[i,1] = unsafe_load(mesh.edgelist, 2*i-1);
            this.edge[i,2] = unsafe_load(mesh.edgelist, 2*i);
            this.edge_marker[i] = unsafe_load(mesh.edgemarkerlist, i);
        end
        this.edge = this.edge + 1;
                                
        return this
    end
end # end type


# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------


function mesh_unit_square(n_segs_per_edge :: Int64, angle = 0.0 :: Float64)

    area_max = round(sqrt(3)/4, 10)/n_segs_per_edge^2
    if angle == 0.0
        switches = string("pzYqa", @sprintf("%10.10f",area_max), "en")
    else
        switches = string("pzYq", angle, "a", @sprintf("%10.10f",area_max), "en")
    end
    println(switches)
    
    mesh_buffer = ccall((:tesselate_unit_square, "libtesselate"), Triangle_mesh_C, (Int64, Cstring), n_segs_per_edge, switches);

    mesh = TriMesh(mesh_buffer)
    mesh_buffer = [];

    return mesh
end


# ---------------------------------------------------------------------------------------------


function mesh_unit_simplex(h = 0.1 :: Float64, angle = 0.0 :: Float64)

    area_max = round(sqrt(3)/4, 10)*h^2
    if angle == 0.0
        switches = string("pzqa", @sprintf("%10.10f",area_max), "en")
    else
        switches = string("pzq", angle, "a", @sprintf("%10.10f",area_max), "en")
    end
    println(switches)
    
    mesh_buffer = ccall((:tesselate_unit_simplex, "libtesselate"), Triangle_mesh_C, (Float64, Cstring), h, switches);

    mesh = TriMesh(mesh_buffer)
    mesh_buffer = [];

    return mesh
end


# ---------------------------------------------------------------------------------------------


function write_mesh(mesh :: TriMesh, filename)

    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh.point)

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[i,:])
        cells = push!(cells, temp)
    end

    vtkfile = vtk_grid(filename, points, cells)

    outfiles = vtk_save(vtkfile)

end


# ---------------------------------------------------------------------------------------------


function plot_mesh(m :: TriMesh)

    fig = matplotlib[:pyplot][:figure]("2D Mesh Plot", figsize = (10,10))

    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")

    # Connectivity list -1 for Python
    tri = ax[:triplot](m.point[:,1], m.point[:,2], m.cell - 1 )
    setp(tri, linestyle = "-",
         marker = "None",
         linewidth = 2,
         color = "green")
    
    fig[:canvas][:draw]()

    return nothing
end


# ---------------------------------------------------------------------------------------------


function get_cell( m :: TriMesh, ind_c = 1:m.n_cell  )
    ind_c = vec(collect(ind_c))
    
    return m.cell[ind_c,:]
end # end get_cell


# ---------------------------------------------------------------------------------------------


function get_point( m :: TriMesh, ind_p = 1:m.n_point )
    
    return m.point[ind_p,:]
end # end get_point


# ---------------------------------------------------------------------------------------------


function map_ref_point(m :: TriMesh, x :: Array{Float64,2}, ind_c = 1:m.n_cell)
    ind_c = vec(collect(ind_c))

    y = transpose( [x  ones(size(x,1))] )
    n_points = size(x,1)
    n_cells = length(ind_c)

    s = "sparse($y[:,:]),"^n_cells
    s1 = s[1:end-1]
    s2 = string("blkdiag(", s1, ")")
    e = parse(s2)
    D = eval(e)

    s = ""
    for i=1:n_cells
        T = m.T_ref2cell[:,:,ind_c[i]]
        s = string(s, "$T,")
    end
    s1 = s[1:end-1]
    s2 = string("hcat(", s1, ")")
    e = parse(s2)
    A = eval(e)

    P0 = full(A*D)
    P0 = reshape(P0, 2, n_points, n_cells)
    P1 = reshape(P0[1,:], n_points, 1, n_cells)
    P2 = reshape(P0[2,:], n_points, 1, n_cells)
    P = hcat(P1, P2)
    
    return P
end


# ---------------------------------------------------------------------------------------------
end # end module Mesh
