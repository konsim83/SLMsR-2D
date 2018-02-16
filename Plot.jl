module Plot

using Mesh, PyPlot

function plot_mesh(m :: TriangleMesh.TriMesh; 
                        linewidth :: Real = 1.5, 
                        marker :: String = "None",
                        markersize :: Real = 5,
                        linestyle :: String = "-",
                        color :: String = "blue")

    fig = matplotlib[:pyplot][:figure]("2D Mesh Plot", figsize = (10,10))
    
    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")
    
    # Connectivity list -1 for Python
    tri = ax[:triplot](m.point[:,1], m.point[:,2], m.cell-1 )
    setp(tri,   linestyle = linestyle,
                linewidth = linewidth,
                marker = marker,
                markersize = markersize,
                color = color)
    
    fig[:canvas][:draw]()
    
    return fig
end



function plot_mesh(m :: TriangleMesh.TriMesh, idx :: Array{Int,1}; 
                        linewidth :: Real = 1.5, 
                        marker :: String = "None",
                        markersize :: Real = 200,
                        linestyle :: String = "-",
                        color :: String = "blue",
                        color_text :: String = "green",
                        offset :: Float64 = 0.02)

    fig = matplotlib[:pyplot][:figure]("2D Mesh Plot", figsize = (10,10))
    
    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")
    
    # Connectivity list -1 for Python
    tri = ax[:triplot](m.point[:,1], m.point[:,2], m.cell-1 )
    setp(tri,   linestyle = linestyle,
                linewidth = linewidth,
                marker = marker,
                color = color)

    scatter(m.point[:,1], m.point[:,2], s=markersize, color="green")

    for i=1:m.n_point
        ax[:text](m.point[i,1]+offset, m.point[i,2]+offset, string(idx[i]), fontsize=20, color=color_text)
    end

    
    fig[:canvas][:draw]()
    
    return fig
end



function plot_mesh_collection(M :: TriMesh_collection, ind_c = 1:M.mesh.n_cell)

    fig = matplotlib[:pyplot][:figure]("2D Mesh collection Plot", figsize = (10,10))
    
    ax = matplotlib[:pyplot][:axes]()
    ax[:set_aspect]("equal")
    
    # Connectivity list -1 for Python
    tri = ax[:triplot](M.mesh.point[:,1], M.mesh.point[:,2], M.mesh.cell - 1 )
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

end