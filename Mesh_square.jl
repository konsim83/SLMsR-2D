function mesh_unit_square(n_segs_per_edge :: Int64, angle = 0.0 :: Float64)

    area_max = round(sqrt(3)/4, 10)/n_segs_per_edge^2
    if angle == 0.0
        switches = string("pzYqa", @sprintf("%10.10f",area_max), "en")
    else
        switches = string("pzYq", angle, "a", @sprintf("%10.10f",area_max), "en")
    end
    #println(switches)
    
    mesh_buffer = ccall((:tesselate_unit_square, "libtesselate"), Triangle_mesh_C, (Int64, Cstring), n_segs_per_edge, switches);

    mesh = TriMesh(mesh_buffer, "unit square")
    mesh_buffer = [];

    return mesh
end
