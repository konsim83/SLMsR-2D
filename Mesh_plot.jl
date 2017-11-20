function plot_mesh(m :: TriMesh)

    fig = matplotlib[:pyplot][:figure]("2D Mesh Plot", figsize = (10,10))
    
    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")
    
    # Connectivity list -1 for Python
    tri = ax[:triplot](m.point[:,1], m.point[:,2], m.cell - 1 )
    setp(tri, linestyle = "-",
         marker = "None",
         linewidth = 2,
         color = "green")
    
    fig[:canvas][:draw]()
    
    return nothing
end



function plot_mesh_collection(M :: TriMesh_collection, ind_c = 1:M.mesh_c.n_cell)

    fig = matplotlib[:pyplot][:figure]("2D Mesh collection Plot", figsize = (10,10))
    
    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")
    
    # Connectivity list -1 for Python
    tri = ax[:triplot](M.mesh_c.point[:,1], M.mesh_c.point[:,2], M.mesh_c.cell - 1 )
    setp(tri, linestyle = "-",
         marker = "None",
         linewidth = 3,
         color = "green")

    for i in ind_c
        # Connectivity list -1 for Python
        tri = ax[:triplot](M.mesh_f[i].point[:,1], M.mesh_f[i].point[:,2], M.mesh_f[i].cell - 1 )
        setp(tri, linestyle = "--",
             marker = "None",
             linewidth = 0.7,
             color = "red")
    end # end for
    
    fig[:canvas][:draw]()
    
    return nothing
    
end