function mesh_triangle_uniform_edges(point:: Array{Float64,2}, n_segs_per_edge :: Int64, angle = 0.0 :: Float64)

    size(point)!=(3,2) ? error("Size of point array must be (3,2).") :
    
    h = minimum([norm(point[1,:]-point[2,:])   norm(point[2,:]-point[3,:])   norm(point[3,:]-point[1,:])]) / n_segs_per_edge
    area_max = round(sqrt(3)/4, 10)*h^2
    
    if angle == 0.0
        switches = string("pzYQqa", @sprintf("%10.10f",area_max), "en")
    else
        switches = string("pzYQq", angle, "a", @sprintf("%10.10f",area_max), "en")
    end
    
    mesh_buffer = ccall((:tesselate_triangle_uniform_edges, "libtesselate"), Triangle_mesh_C,
                        (Float64, Float64, Float64, Float64, Float64, Float64, Int64, Cstring),
                        point[1,1], point[1,2], point[2,1], point[2,2], point[3,1], point[3,2], n_segs_per_edge, switches);

    mesh = TriMesh(mesh_buffer, "Triangle - number of points per boundary edge is fixed")
    mesh_buffer = [];

    return mesh
end
