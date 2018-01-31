function mesh_unit_simplex(n_refinements :: Int64)

    point = [0.0 0.0 ; 
            1.0 0.0 ; 
            0.0 1.0]

    mesh = mesh_triangle(point, 1)

    for i=1:n_refinements
        mesh = TriangleMesh.refine_rg(mesh)
    end

    return mesh
end


function mesh_unit_simplex(n_segs_per_edge :: Int64)
    
    point = [0.0 0.0 ; 
            1.0 0.0 ; 
            0.0 1.0]

    return mesh_triangle(point, n_segs_per_edge)
end


# ---------------------------------------------------------------------------------------------------


function mesh_triangle(point:: Array{Float64,2}, n_segs_per_edge :: Int64)

    #########################################
    #########################################

    return mesh
end


# ---------------------------------------------------------------------------------------------------


function mesh_unit_square(n_segs_per_edge :: Int64)

    #########################################
    #########################################

    return mesh
end


# ---------------------------------------------------------------------------------------------------



function refine_mesh_unit_square(mesh :: TriMesh, n_refinements :: Int64)

    for i=1:n_refinements
        mesh = TriangleMesh.refine_rg(mesh)
    end

    return mesh
end


# ---------------------------------------------------------------------------------------------------