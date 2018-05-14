"""
    struct Solution_FEM <: AbstractSolution

    Solution type for standard FEM.

"""
struct Solution_FEM <: AbstractSolution
    # Solution at nodes
    u :: Array{Float64,2}

    function Solution_FEM(u :: Array{Float64,2})

        return new(u)
    end
end # end type


# --------------------------------------------------------------------
# --------------------------------------------------------------------


# Outer constructor
function Solution_FEM(dof :: AbstractDof,
                          par :: Parameter.Parameter_FEM)

    # Reserve memory for the solution
    u = Array{Float64,2}(dof.n_node, par.n_steps+1)
    
    return Solution_FEM(u)
end


"""
    Solution_FEM(u_in :: Array{Array{Array{Float64,1},1},1},
                        n_node :: Int64)

    Outer constructor for struct 'Solution_FEM' from post-processed FEM data
    (solution evaluated at a suitable set of points).

"""
function Solution_FEM(u_in :: Array{Array{Array{Float64,1},1},1},
                        n_node :: Int64)

    # Reserve memory for the solution
    u = Array{Float64,2}(n_node, length(u_in[1][1][:]))

    for i in 1:length(u_in)
        u[i,:] = u_in[i][1][:]
    end

    return Solution_FEM(u)
end


"""
    Solution_FEM(u_in :: Array{Array{Array{Array{Float64,1},1},1},1},
                        n_node :: Int64)

    Outer constructor for struct 'Solution_FEM' from post-processed MsFEM data
    (solution evaluated at a suitable set of points).

"""
function Solution_FEM(u_in :: Array{Array{Array{Array{Float64,1},1},1},1},
                        n_node :: Int64)

    # Reserve memory for the solution
    u = Array{Float64,2}(n_node, length(u_in[1][1][1][:]))

    for i in 1:length(u_in)
        u[i,:] = u_in[i][1][1][:]
    end

    return Solution_FEM(u)
end