function refine_rg(mesh :: TriangleMesh.TriMesh, n_refinements :: Int64)

    for i=1:n_refinements
        mesh = TriangleMesh.refine_rg(mesh)
    end

    return mesh
end

# ---------------------------------------------------------------------------------------------------