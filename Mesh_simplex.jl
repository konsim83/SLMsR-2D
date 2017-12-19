function mesh_unit_simplex(h = 0.1 :: Float64, angle = 33.0 :: Float64)

    area_max = round(sqrt(3)/4, 10)*h^2
    if angle == 0.0
        switches = string("pzYqa", @sprintf("%10.10f",area_max), "en")
    else
        switches = string("pzYq", angle, "a", @sprintf("%10.10f",area_max), "en")
    end
    
    mesh_buffer = ccall((:tesselate_unit_simplex, "libtesselate"), Triangle_mesh_C, (Float64, Cstring), h, switches);

    mesh = TriMesh(mesh_buffer, "Simplex")
    mesh_buffer = [];

    return mesh
end


# ---------------------------------------------------------------------------------------------------


function mesh_unit_simplex_uniform_edges(n_segs_per_edge :: Int64, angle = 0.0 :: Float64)

    h = 1/n_segs_per_edge
    area_max = round(sqrt(3)/4, 10)*h^2
    
    if angle == 0.0
        switches = string("pzYqa", @sprintf("%10.10f",area_max), "en")
    else
        switches = string("pzYq", angle, "a", @sprintf("%10.10f",area_max), "en")
    end
    
    mesh_buffer = ccall((:tesselate_unit_simplex_uniform_edges, "libtesselate"), Triangle_mesh_C, (Int64, Cstring), n_segs_per_edge, switches);

    mesh = TriMesh(mesh_buffer, "Simplex - number of points per boundary edge is fixed")
    mesh_buffer = [];

    return mesh
end
