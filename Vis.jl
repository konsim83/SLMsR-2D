module Vis

using FEM, Mesh, WriteVTK

function writeSolution(solution :: FEM.Solution_FEM, mesh :: Mesh.TriMesh, filename)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh.point)

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[i,:])
        cells = push!(cells, temp)
    end
    
    vtkfile = vtk_grid(filename, points, cells)

    vtk_point_data(vtkfile, convert(Array{Float32}, solution.u[:,1]), "point data")
    
    outfiles = vtk_save(vtkfile)

    return nothing
end

end # end module
