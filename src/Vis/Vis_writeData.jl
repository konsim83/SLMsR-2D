function writeVelocityField(problem :: Problem.AbstractProblem,
							mesh :: Mesh.TriangleMesh.TriMesh,
							T :: Array{Float64,1},
							filename :: String)

	# Write files into #PWD/meshfiles folder
	if ~ispath(pwd() * "/data/bgVel")
		mkdir("data/bgVel")
	end
	path = pwd() * "/data/bgVel/"  
	filename = path * basename(filename)

	cell_type = VTKCellTypes.VTK_TRIANGLE
    points = mesh.point

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[:,i])
        cells = push!(cells, temp)
    end

    N = length(T)
    p = Progress(N, 0.01, "Writing vecolicty field ...", 10)
    
    k_time = 1
    for t in T
    	velocity = Problem.velocity(problem, t, mesh.point)
        vtkfile = vtk_grid(string(filename, "_", lpad(k_time,4,"0")), points, cells)
        # vtk_point_data(vtkfile, convert(Array{Float32}, solution.u[:,k_time]), "point data")
        vtk_point_data(vtkfile, hcat(velocity...), "bg_vel")
    
        outfiles = vtk_save(vtkfile)

        k_time += 1
        next!(p)
    end
end