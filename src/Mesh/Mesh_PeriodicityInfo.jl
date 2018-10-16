struct PeriodicityInfo

	edge_marker_pair :: Array{Int,2}

	# Contains info on how to map the second edge to the first for each row in
	# 'edge_marker_pair'.
	f_edge_to_edge :: Array{Function,1}

end

function DoublePeriodicUnitSquare()

	# Create the maps that translate dofs to mesh based
    # variables. This is one of the core modules for the
    # handling DoFs of periodic meshes (without actually
    # touching the mesh).

    # Tells where any mesh based variable can be found
    # topologically

    # Needs information on the mesh
    edge_marker_pair = [1 3 ; 2 4]
    f_edge_to_edge = Array{Function,1}(undef, 2)
    f_edge_to_edge[1] = function(p :: Array{Float64,2})
        # Periodicity in y-direction, maps upper edge to lower edge

        return broadcast(+, [0.0 ;-1.0], p)
    end

    f_edge_to_edge[2] = function(p :: Array{Float64,2})
        # Periodicity in x-direction, maps left edge to right edge

        return broadcast(+, [1.0 ;0.0], p)
    end

    return PeriodicityInfo(edge_marker_pair, f_edge_to_edge)
end


function X_PeriodicUnitSquare()

    # Create the maps that translate dofs to mesh based
    # variables. This is one of the core modules for the
    # handling DoFs of periodic meshes (without actually
    # touching the mesh).

    # Tells where any mesh based variable can be found
    # topologically

    # Needs information on the mesh
    edge_marker_pair = [2 4]
    f_edge_to_edge = Array{Function,1}(1)

    f_edge_to_edge[1] = function(p :: Array{Float64,2})
        # Periodicity in x-direction, maps left edge to right edge

        return broadcast(+, [1.0 ;0.0], p)
    end

    return PeriodicityInfo(edge_marker_pair, f_edge_to_edge)
end


function Y_PeriodicUnitSquare()

    # Create the maps that translate dofs to mesh based
    # variables. This is one of the core modules for the
    # handling DoFs of periodic meshes (without actually
    # touching the mesh).

    # Tells where any mesh based variable can be found
    # topologically

    # Needs information on the mesh
    edge_marker_pair = [1 3]
    f_edge_to_edge = Array{Function,1}(1)
    f_edge_to_edge[1] = function(p :: Array{Float64,2})
        # Periodicity in y-direction, maps upper edge to lower edge

        return broadcast(+, [0.0 ;-1.0], p)
    end

    return PeriodicityInfo(edge_marker_pair, f_edge_to_edge)
end