include("Mesh.jl")


# ----------------------------------------
function identify_points(mesh :: Mesh.TriangleMesh.TriMesh,
                            edge_marker_pair :: Array{Int64,2},
                            edge_trafo :: Array{Function,1})

    length(edge_trafo)!=size(edge_marker_pair,1) ? 
        error("Number of Matched edge pairs must coincide with transformations.") :

    # index_map = collect(1:mesh.n_point)
    ind_point_boundary = find(mesh.point_marker.!=0)
    index_map = []
    for i in ind_point_boundary
        push!(index_map,[i])
    end

    for i in 1:size(edge_marker_pair,1)
        pair = edge_marker_pair[i,:]

        point_ind_on_edge1 = sort(unique(vec(mesh.edge[mesh.edge_marker.==pair[1],:]')))
        # display(point_ind_on_edge1)
        point_ind_on_edge2 = sort(unique(vec(mesh.edge[mesh.edge_marker.==pair[2],:]')))

        f = edge_trafo[i]

        point_edge1 = mesh.point[point_ind_on_edge1,:]
        # Apply f to edge1, gives permuted edge2
        point_edge2_on_edge1 = f(mesh.point[point_ind_on_edge2,:])

        # Find the indices of mapped points on edge2 in edge1, then change the
        # index of points on edge2 to get identified with edge1
        for j=1:size(point_edge2_on_edge1,1)
            ind_p2_in_p1 = closest_index(point_edge1, point_edge2_on_edge1[j,:])
            push!(index_map[point_ind_on_edge2[ind_p2_in_p1]],point_ind_on_edge1[j])
        end
    end
    # index_map = map(x->unique(sort(x)), index_map)

    # for k=1:length(index_map)
    #     while length(index_map[k])>1
    #         ind_2b_replaced = index_map[k][end]
    #         ind_new = index_map[k][1]
            
    #         map(x-> (x[x.==ind_2b_replaced]=ind_new), index_map);
           
    #         # Keep only unique values
    #         index_map = map(x->unique(x), index_map)
    #     end        
    # end

    # index_map = vcat(index_map...)

    # # Squeze vector
    # d = setdiff(1:maximum(index_map), index_map)
    # while length(d)>0
    #     index_map[index_map.>d[1]] -= 1
    #     d = setdiff(1:maximum(index_map), index_map)
    # end
    
    # # append the inner points
    # n_point_interior = sum(mesh.point_marker.==0)
    # ind_point_interior = collect(1:n_point_interior) + maximum(index_map)
    # append!(index_map, ind_point_interior)

    return index_map
end

# find closest point in array
function closest_index(P :: Array{Float64,2}, p :: Array{Float64,1})

    ibest = start(eachindex(P[:,1]))
    dxbest = sum(abs.(P[ibest,:] - p))

    for ind in eachindex(P[:,1])
        dx = sum(abs.(P[ind,:] - p))
        if dx < dxbest
            dxbest = dx
            ibest = ind
        end
    end
    
    return ibest
end
# ----------------------------------------




mesh = Mesh.mesh_unit_square(n)
mesh = Mesh.refine_rg(mesh, 1)

edge_marker_pair = [1 3 ; 2 4]
f_edge_to_edge = Array{Function,1}(2)
f_edge_to_edge[1] = function(p :: Array{Float64,2})
    # Periodicity in y-direction, maps upper edge to lower edge

    return broadcast(+, [0.0 -1.0], p)
end

f_edge_to_edge[2] = function(p :: Array{Float64,2})
    # Periodicity in x-direction, maps left edge to right edge

    return broadcast(+, [1.0 0.0], p)
end
@time index_map = identify_points(mesh, edge_marker_pair, f_edge_to_edge)