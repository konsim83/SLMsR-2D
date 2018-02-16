# Actual Julia mesh type for collection of meshes (usefull for
# multiscale and multigrid)
mutable struct TriMesh_collection

    mesh :: TriangleMesh.TriMesh
    mesh_f :: Array{TriangleMesh.TriMesh, 1}

    n_elem_f :: Array{Int64, 1}
    
    function TriMesh_collection(mesh :: TriangleMesh.TriMesh, mesh_simplex :: TriangleMesh.TriMesh)
        this = new()

        this.mesh = mesh
        this.mesh_f = Array{TriangleMesh.TriMesh}(mesh.n_cell)
        this.n_elem_f = Array{Int64}(mesh.n_cell)

        for i=1:mesh.n_cell
            T_ref2cell = [mesh.point[mesh.cell[i,2],:]-mesh.point[mesh.cell[i,1],:] mesh.point[mesh.cell[i,3],:]-mesh.point[mesh.cell[i,1],:]  mesh.point[mesh.cell[i,1],:]]

            point_mapped = [mesh_simplex.point ones(size(mesh_simplex.point,1))] * T_ref2cell'

            this.mesh_f[i] = TriangleMesh.TriMesh(mesh_simplex, point_mapped, string("Mapped ", mesh_simplex.mesh_info))
            this.n_elem_f[i] = this.mesh_f[i].n_cell
        end

        return this
    end # end constructor

    function TriMesh_collection(mesh :: TriangleMesh.TriMesh, n_segs_per_edge_f :: Int64)
        this = new()

        this.mesh = mesh
        this.mesh_f = Array{TriangleMesh.TriMesh}(mesh.n_cell)
        this.n_elem_f = Array{Int64}(mesh.n_cell)

        point = get_point(mesh, get_cell(mesh, 1:mesh.n_cell))
        for i=1:mesh.n_cell
            this.mesh_f[i] = mesh_triangle(point[:,:,i], n_segs_per_edge_f)
            this.n_elem_f[i] = this.mesh_f[i].n_cell
        end

        return this
    end # end constructor

end # end type