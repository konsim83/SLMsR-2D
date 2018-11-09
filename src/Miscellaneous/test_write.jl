using ProgressMeter, WriteVTK, Mesh

mesh = Mesh.mesh_unit_square(100)
t = collect(range(0, stop=1, length=101))

function writeVectorField(field :: Function,
							mesh :: Mesh.TriangleMesh.TriMesh,
							TT :: Array{Float64,1},
							filename = "field" :: String)

	# Write files into #PWD/meshfiles folder
	if ~ispath(pwd() * "/data/testData")
		mkdir("data/testData")
	end
	path = pwd() * "/data/testData/"
	filename = path * basename(filename)

	cell_type = VTKCellTypes.VTK_TRIANGLE
    points = mesh.point

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[:,i])
        cells = push!(cells, temp)
    end

    N = length(TT)
    p = Progress(N, 0.01, "Writing test field ...", 10)
    
    k_time = 1
    for tt in TT
    	velocity = field(tt, mesh.point)
        vtkfile = vtk_grid(string(filename, "_", lpad(k_time,4,"0")), points, cells)
        # vtk_point_data(vtkfile, convert(Array{Float32}, solution.u[:,k_time]), "point data")
        vtk_point_data(vtkfile, hcat(velocity...), "TestField")
    
        outfiles = vtk_save(vtkfile)

        k_time += 1
        next!(p)
    end
end


function writeScalarField(scalar :: Function,
                            mesh :: Mesh.TriangleMesh.TriMesh,
                            TT :: Array{Float64,1},
                            filename = "scalar" :: String)

    # Write files into #PWD/meshfiles folder
    if ~ispath(pwd() * "/data/testData")
        mkdir("data/testData")
    end
    path = pwd() * "/data/testData/"
    filename = path * basename(filename)

    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = mesh.point

    cells = MeshCell[]
    for i=1:mesh.n_cell
        temp = MeshCell(cell_type, mesh.cell[:,i])
        cells = push!(cells, temp)
    end

    N = length(TT)
    p = Progress(N, 0.01, "Writing test scalar ...", 10)
    
    k_time = 1
    for tt in TT
        u = scalar(tt, mesh.point)
        vtkfile = vtk_grid(string(filename, "_", lpad(k_time,4,"0")), points, cells)
        # vtk_point_data(vtkfile, convert(Array{Float32}, solution.u[:,k_time]), "point data")
        vtk_point_data(vtkfile, u, "TestScalar")
    
        outfiles = vtk_save(vtkfile)

        k_time += 1
        next!(p)
    end
end


function scalar(tt :: Float64, x :: Array{Float64,2})
    k1 = 1
    k2 = 1

    return sin.(2*k1*pi*(x[1,:].-tt)) .* cos.(2*k2*pi*x[2,:])
end 


function field(tt :: Float64, x :: Array{Float64,2})
    k1 = 1
    k2 = 1

    V = hcat(-sin.(2*pi*k1*(x[1,:].-tt)) .* sin.(2*pi*k2*(x[2,:])) *2*pi*k2,
                -cos.(2*pi*k1*(x[1,:].-tt)) .* cos.(2*pi*k2*(x[2,:])) *2*pi*k1
            )
    return [ V[i,:] for i in 1:size(x,2) ]
end 