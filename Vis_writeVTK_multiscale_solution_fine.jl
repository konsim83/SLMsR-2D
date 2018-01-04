function writeSolution_at_timestep(solution :: FEM.Solution_MsFEM, mesh_collection :: Mesh.TriMesh_collection, k_time :: Int64, filename :: String)

    P, C = unify_meshes(mesh_collection, solution)
    D = unify_data(mesh_collection, solution, k_time)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(P)

    cells = MeshCell[]
    for i=1:size(C,1)
        temp = MeshCell(cell_type, C[i,:])
        cells = push!(cells, temp)
    end
    
    vtkfile = vtk_grid(string("./data/", filename, "_", lpad(k_time,4,0)), points, cells)

    vtk_point_data(vtkfile, convert(Array{Float32}, D), "point data")
    
    outfiles = vtk_save(vtkfile)

    return nothing
end


# ================================================


function writeSolution_all(solution :: FEM.Solution_MsFEM, mesh_collection :: Mesh.TriMesh_collection, filename :: String)
    
    P, C = unify_meshes(mesh_collection, solution)
    
    cell_type = VTKCellTypes.VTK_TRIANGLE
    points = transpose(P)

    cells = MeshCell[]
    for i=1:size(C,1)
        temp = MeshCell(cell_type, C[i,:])
        cells = push!(cells, temp)
    end

    for k_time in 1:size(solution.u,2)
        vtkfile = vtk_grid(string("./data/", filename, "_", lpad(k_time,4,0)), points, cells)
        vtk_point_data(vtkfile, convert(Array{Float32}, unify_data(mesh_collection, solution, k_time)), "point data")
    
        outfiles = vtk_save(vtkfile)
    end

    return nothing
end


# ================================================


# -------   Outer constructor for TriMesh   -------
function unify_meshes(mesh_collection :: Mesh.TriMesh_collection, solution :: FEM.Solution_MsFEM)

    P = mesh_collection.mesh_f[1].point
    C = mesh_collection.mesh_f[1].cell
    
    for i in 2:length(mesh_collection.mesh_f)
        P, C = add_meshes(P, C, mesh_collection.mesh_f[i])
    end

    return P, C
end


function unify_data(mesh_collection :: Mesh.TriMesh_collection, solution :: FEM.Solution_MsFEM, k_time)

    D = gather_data(mesh_collection, solution, 1, k_time)
    
    for i in 2:length(mesh_collection.mesh_f)
        D = vcat(D, gather_data(mesh_collection, solution, i, k_time))
    end

    return D
end


function gather_data(mesh_collection :: Mesh.TriMesh_collection, solution :: FEM.Solution_MsFEM, ind_cell, k_time)

    mesh = mesh_collection.mesh_f[ind_cell]
    cell_coarse = mesh_collection.mesh.cell[ind_cell,:]
    
    phi1 = solution.phi_1[ind_cell][:,k_time]
    phi2 = solution.phi_2[ind_cell][:,k_time]
    phi3 = solution.phi_3[ind_cell][:,k_time]

    u = solution.u[cell_coarse,k_time]

    data = u[1]*phi1 + u[2]*phi2 + u[3]*phi3
    
    return data
end

# -------   Add two meshes   -------
function add_meshes(p1 :: Array{Float64,2}, c1 :: Array{Int64,2}, mesh2 :: Mesh.TriMesh)
    
    P = vcat(p1, mesh2.point)
    C = vcat(c1, mesh2.cell+size(p1,1))
    
    #=
    p2, c2 = copy(mesh2.point), copy(mesh2.cell)
    
    n_p1 = size(p1,1)
    n_p2 = mesh2.n_point
    
    n_found = 0
    ind_delete = []
    for j in 1:n_p1
        #p1_in_p2 = (Bool[ p1[j,:] == p2[i,:] for i in 1:n_p2 ])
        p1_in_p2 = (Bool[ sum(abs.(p1[j,:] - p2[i,:]))<1e-10 for i in 1:n_p2 ])
        
        # This index in c2 must be set to j
        ind_p1_in_p2 = find(p1_in_p2)
        if ~isempty(ind_p1_in_p2)
            push!(ind_delete, ind_p1_in_p2[1])
            c2[c2.==ind_p1_in_p2] = -j # set it to -j to easily identify changed points
            # increase by one if a point is found
            n_found += sum(p1_in_p2)
        end
    end
    
    p2 = p2[setdiff(1:n_p2,ind_delete),:]
    c2[c2.>0] = c2[c2.>0]  + n_p1 - n_found 
    c2 = abs.(c2)
    
    P = vcat(p1, p2)
    C = vcat(c1, c2)
    
    C = unique(C,1)
    =#
    
    return P, C
end
