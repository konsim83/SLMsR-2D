function writeSolution_at_timestep(solution :: FEM.Solution_FEM, mesh :: Mesh.TriMesh, k_time :: Int64, filename :: String)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh.point)

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[i,:])
        cells = push!(cells, temp)
    end
    
    vtkfile = vtk_grid(string("./data/", filename, "_", lpad(k_time,4,0)), points, cells)

    vtk_point_data(vtkfile, convert(Array{Float32}, solution.u[:,k_time]), "point data")
    
    outfiles = vtk_save(vtkfile)

    return nothing
end


# ================================================


function writeSolution_all(solution :: FEM.Solution_FEM, mesh :: Mesh.TriMesh, filename :: String)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh.point)

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[i,:])
        cells = push!(cells, temp)
    end

    for k_time in 1:size(solution.u,2)
        vtkfile = vtk_grid(string("./data/", filename, "_", lpad(k_time,4,0)), points, cells)
        vtk_point_data(vtkfile, convert(Array{Float32}, solution.u[:,k_time]), "point data")
    
        outfiles = vtk_save(vtkfile)
    end

    return nothing
end
