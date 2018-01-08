function writeBasis_at_timestep(solution :: FEM.Solution_MsFEM,
                                mesh_collection :: Mesh.TriMesh_collection,
                                k_time :: Int64,
                                ind_cell :: Int64,
                                filename :: String)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh_collection.mesh_f[ind_cell].point)

    cells = MeshCell[]
    for i=1:mesh_collection.mesh_f[ind_cell].n_cell
        temp = MeshCell(cell_type, mesh_collection.mesh_f[ind_cell].cell[i,:])
        cells = push!(cells, temp)
    end
    
    vtkfile = vtk_grid(string("./data/", filename, "_cell-", string(ind_cell), "_", lpad(k_time,4,0)), points, cells)

    vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_1[ind_cell][:,k_time]), "basis 1 - nodal values")
    vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_2[ind_cell][:,k_time]), "basis 2 - nodal values")
    vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_3[ind_cell][:,k_time]), "basis 3 - nodal values")
    
    outfiles = vtk_save(vtkfile)

    return nothing
end


# ==============================================================


function writeBasis_all_steps(solution :: FEM.Solution_MsFEM,
                              mesh_collection :: Mesh.TriMesh_collection,
                              ind_cell :: Int64,
                              filename :: String)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(mesh_collection.mesh_f[ind_cell].point)

    cells = MeshCell[]
    for i=1:mesh_collection.mesh_f[ind_cell].n_cell
        temp = MeshCell(cell_type, mesh_collection.mesh_f[ind_cell].cell[i,:])
        cells = push!(cells, temp)
    end

    N = size(solution.u,2)
    p = Progress(N, 0.01, "Writing progress of one basis function ...", 10)

    for k_time in 1:size(solution.u,2)
        vtkfile = vtk_grid(string("./data/", filename, "_cell-", string(ind_cell), "_", lpad(k_time,4,0)), points, cells)
        
        vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_1[ind_cell][:,k_time]), "basis 1 - nodal values")
        vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_2[ind_cell][:,k_time]), "basis 2 - nodal values")
        vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_3[ind_cell][:,k_time]), "basis 3 - nodal values")
    
        outfiles = vtk_save(vtkfile)

        next!(p)
    end

    return nothing
end


# ==============================================================


function writeBasis_all(solution :: FEM.Solution_MsFEM,
                        mesh_collection :: Mesh.TriMesh_collection,
                        filename :: String)

    N = mesh_collection.mesh.n_cell
    p = Progress(N, 0.01, "Writing progress of all basis functions ...", 10)
    
    for ind_cell = 1:mesh_collection.mesh.n_cell
        cell_type = VTKCellTypes.VTK_TRIANGLE
        points = transpose(mesh_collection.mesh_f[ind_cell].point)

        cells = MeshCell[]
        for i=1:mesh_collection.mesh_f[ind_cell].n_cell
            temp = MeshCell(cell_type, mesh_collection.mesh_f[ind_cell].cell[i,:])
            cells = push!(cells, temp)
        end
        
        for k_time in 1:size(solution.u,2)
            vtkfile = vtk_grid(string("./data/", filename, "_cell-", string(ind_cell), "_", lpad(k_time,4,0)), points, cells)
            
            vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_1[ind_cell][:,k_time]), "basis 1 - nodal values")
            vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_2[ind_cell][:,k_time]), "basis 2 - nodal values")
            vtk_point_data(vtkfile, convert(Array{Float32}, solution.phi_3[ind_cell][:,k_time]), "basis 3 - nodal values")
            
            outfiles = vtk_save(vtkfile)
        end

        next!(p)
    end

    return nothing
end
