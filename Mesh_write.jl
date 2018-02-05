function write_mesh_vtk(mesh :: TriangleMesh.TriMesh, filename)

    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh.point)
    
    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[i,:])
        cells = push!(cells, temp)
    end
    
    vtkfile = vtk_grid(filename, points, cells)
            
    outfiles = vtk_save(vtkfile)

    return nothing
end
