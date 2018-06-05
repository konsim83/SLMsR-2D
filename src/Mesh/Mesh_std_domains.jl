# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------

function mesh_triangle(point:: Array{Float64,2}, n_segs_per_edge :: Int64)

    n_segment = n_segs_per_edge * 3
    n_point = n_segs_per_edge * 3
    n_point_marker = 1
    n_point_attribute = 0
    n_hole = 0

    p = TriangleMesh.Polygon_pslg(n_point, n_point_marker, n_point_attribute, n_segment, n_hole)

    segment = zeros(Int64, n_segment,2)
    segment_marker = zeros(Int64, n_segment)
    point = zeros(n_point,2)
    point_marker = zeros(Int64,n_point,1)

    # Add the points, edge 1
    segment[end,:] = [n_point ; 1]
    segment_marker[end] = 3
    point_marker[1,1] = 1
    point[1,:] = [0.0 ; 0.0]
    for i = 2:(n_segs_per_edge)
        point[i,:] = point[i-1,:] + [1.0 ; 0.0]/n_segs_per_edge
        point_marker[i,1] = 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 1
    end

    # Add the points, edge 2
    segment[n_segs_per_edge,:] = [n_segs_per_edge ; n_segs_per_edge+1]
    segment_marker[n_segs_per_edge] = 1
    point_marker[n_segs_per_edge+1,1] = 1
    point[n_segs_per_edge+1,:] = [1.0 ; 0.0]
    for i in (n_segs_per_edge+2):(2*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [-1.0 ; 1.0]/n_segs_per_edge
        point_marker[i,1] = 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 2
    end

    # Add the points, edge 3
    segment[2*n_segs_per_edge,:] = [2*n_segs_per_edge ; 2*n_segs_per_edge+1]
    segment_marker[2*n_segs_per_edge] = 2
    point_marker[2*n_segs_per_edge+1,1] = 1
    point[2*n_segs_per_edge+1,:] = [0.0 ; 1.0]
    for i in (2*n_segs_per_edge+2):(3*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [0.0 ; -1.0]/n_segs_per_edge
        point_marker[i,1] = 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 3
    end

    TriangleMesh.set_polygon_point!(p, point)
    TriangleMesh.set_polygon_point_marker!(p, point_marker)
    TriangleMesh.set_polygon_segment!(p, segment)
    TriangleMesh.set_polygon_segment_marker!(p, segment_marker)

    switches = "pQYDqen"

    mesh = TriangleMesh.create_mesh(p, switches, info_str = "Delaunay mesh of triangle")

    return mesh
end


function mesh_triangle(point:: Array{Float64,2}, n_segs_per_edge :: Int64, areaMax :: Float64)

    n_segment = n_segs_per_edge * 3
    n_point = n_segs_per_edge * 3
    n_point_marker = 1
    n_point_attribute = 0
    n_hole = 0

    p = TriangleMesh.Polygon_pslg(n_point, n_point_marker, n_point_attribute, n_segment, n_hole)

    segment = zeros(Int64, n_segment,2)
    segment_marker = zeros(Int64, n_segment)
    point = zeros(n_point,2)
    point_marker = zeros(Int64,n_point,1)

    # Add the points, edge 1
    segment[end,:] = [n_point ; 1]
    segment_marker[end] = 3
    point_marker[1,1] = 1
    point[1,:] = [0.0 ; 0.0]
    for i = 2:(n_segs_per_edge)
        point[i,:] = point[i-1,:] + [1.0 ; 0.0]/n_segs_per_edge
        point_marker[i,1] = 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 1
    end

    # Add the points, edge 2
    segment[n_segs_per_edge,:] = [n_segs_per_edge ; n_segs_per_edge+1]
    segment_marker[n_segs_per_edge] = 1
    point_marker[n_segs_per_edge+1,1] = 1
    point[n_segs_per_edge+1,:] = [1.0 ; 0.0]
    for i in (n_segs_per_edge+2):(2*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [-1.0 ; 1.0]/n_segs_per_edge
        point_marker[i,1] = 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 2
    end

    # Add the points, edge 3
    segment[2*n_segs_per_edge,:] = [2*n_segs_per_edge ; 2*n_segs_per_edge+1]
    segment_marker[2*n_segs_per_edge] = 2
    point_marker[2*n_segs_per_edge+1,1] = 1
    point[2*n_segs_per_edge+1,:] = [0.0 ; 1.0]
    for i in (2*n_segs_per_edge+2):(3*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [0.0 ; -1.0]/n_segs_per_edge
        point_marker[i,1] = 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 3
    end

    TriangleMesh.set_polygon_point!(p, point)
    TriangleMesh.set_polygon_point_marker!(p, point_marker)
    TriangleMesh.set_polygon_segment!(p, segment)
    TriangleMesh.set_polygon_segment_marker!(p, segment_marker)

    switches = "pQYDqena" * string(areaMax)

    mesh = TriangleMesh.create_mesh(p, switches, info_str = "Delaunay mesh of triangle")

    return mesh
end


# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------


function mesh_unit_simplex()

    point = [0.0 0.0 ; 
            1.0 0.0 ; 
            0.0 1.0]

    mesh = mesh_triangle(point, 1)

    return mesh
end


function mesh_unit_simplex(n_segs_per_edge :: Int64)
    
    point = [0.0 0.0 ; 
            1.0 0.0 ; 
            0.0 1.0]

    return mesh_triangle(point, n_segs_per_edge)
end


function mesh_unit_simplex(n_segs_per_edge :: Int64, areaMax :: Float64)
    
    point = [0.0 0.0 ; 
            1.0 0.0 ; 
            0.0 1.0]

    return mesh_triangle(point, n_segs_per_edge, areaMax)
end


# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------


function mesh_unit_square()
    
    return mesh_unit_square(1)
end


function mesh_unit_square(n_segs_per_edge :: Int64)

    n_segment = n_segs_per_edge * 4
    n_point = n_segs_per_edge * 4
    n_point_marker = 1
    n_point_attribute = 0
    n_hole = 0

    p = TriangleMesh.Polygon_pslg(n_point, n_point_marker, n_point_attribute, n_segment, n_hole)

    segment = zeros(Int64, n_segment,2)
    segment_marker = zeros(Int64, n_segment)
    point = zeros(n_point,2)
    point_marker = zeros(Int64,n_point,1)

    # Add the points, edge 1
    segment[end,:] = [n_point ; 1]
    segment_marker[end] = 4
    point_marker[1,1] = 1
    point[1,:] = [0.0 ; 0.0]
    for i = 2:(n_segs_per_edge)
        point[i,:] = point[i-1,:] + [1.0 ; 0.0]/n_segs_per_edge
        point_marker[i,1] = 1 #i
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 1
    end

    # Add the points, edge 2
    segment[n_segs_per_edge,:] = [n_segs_per_edge ; n_segs_per_edge+1]
    segment_marker[n_segs_per_edge] = 1
    point_marker[n_segs_per_edge+1,1] = 1
    point[n_segs_per_edge+1,:] = [1.0 ; 0.0]
    for i in (n_segs_per_edge+2):(2*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [0.0 ; 1.0]/n_segs_per_edge
        point_marker[i,1] = 1 #i - 1
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 2
    end

    # Add the points, edge 3
    segment[2*n_segs_per_edge,:] = [2*n_segs_per_edge ; 2*n_segs_per_edge+1]
    segment_marker[2*n_segs_per_edge] = 2
    point_marker[2*n_segs_per_edge+1,1] = 1
    point[2*n_segs_per_edge+1,:] = [1.0 ; 1.0]
    for i in (2*n_segs_per_edge+2):(3*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [-1.0 ; 0.0]/n_segs_per_edge
        point_marker[i,1] = 1 #i
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 3
    end

    # Add the points, edge 4
    segment[3*n_segs_per_edge,:] = [3*n_segs_per_edge ; 3*n_segs_per_edge+1]
    segment_marker[3*n_segs_per_edge] = 3
    point_marker[3*n_segs_per_edge+1,1] = 1
    point[3*n_segs_per_edge+1,:] = [0.0 ; 1.0]
    for i in (3*n_segs_per_edge+2):(4*n_segs_per_edge)
        point[i,:] = point[i-1,:] + [0.0 ; -1.0]/n_segs_per_edge
        point_marker[i,1] = 1 #i
        segment[i-1,:] = [i-1 ; i]
        segment_marker[i-1] = 4
    end

    # Change point markers on edges 3 and 4
    point_marker[(2*n_segs_per_edge+2):(3*n_segs_per_edge)] = reverse(point_marker[2:(n_segs_per_edge)])
    point_marker[(3*n_segs_per_edge+2):(4*n_segs_per_edge)] = reverse(point_marker[(n_segs_per_edge+2):(2*n_segs_per_edge)])

    TriangleMesh.set_polygon_point!(p, point)
    TriangleMesh.set_polygon_point_marker!(p, point_marker)
    TriangleMesh.set_polygon_segment!(p, segment)
    TriangleMesh.set_polygon_segment_marker!(p, segment_marker)

    area_max = 1/n_segs_per_edge^2 * sqrt(3)/4
    switches = "pQYqena" * @sprintf("%.12f", area_max)

    mesh = TriangleMesh.create_mesh(p, switches, info_str = "Delaunay mesh of triangle")

    return mesh
end
# ---------------------------------------------------------------------------------------------------