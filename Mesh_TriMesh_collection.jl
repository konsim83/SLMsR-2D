# Actual Julia mesh type for collection of meshes (usefull for
# multiscale and multigrid)
type TriMesh_collection

    mesh :: TriMesh
    mesh_f :: Array{TriMesh, 1}
    
    function TriMesh_collection(mesh :: TriMesh, mesh_simplex :: TriMesh)
        this = new()

        this.mesh = mesh
        this.mesh_f = Array{TriMesh}(mesh.n_cell)

        point_mapped = map_ref_point(mesh, mesh_simplex.point, 1:mesh.n_cell)
        for i=1:mesh.n_cell
            this.mesh_f[i] = TriMesh(mesh_simplex, point_mapped[:,:,i], string("Mapped ", mesh_simplex.mesh_info))
        end

        return this
    end # end constructor

    function TriMesh_collection(mesh :: TriMesh, n_segs_per_edge_f :: Int64)
        this = new()

        this.mesh = mesh
        this.mesh_f = Array{TriMesh}(mesh.n_cell)

        point = get_point(mesh, get_cell(mesh, 1:mesh.n_cell))
        for i=1:mesh.n_cell
            this.mesh_f[i] = mesh_triangle_uniform_edges(point[:,:,i], n_segs_per_edge_f)
        end

        return this
    end # end constructor

end # end type
