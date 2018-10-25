"""
    struct Solution_FEM <: AbstractSolution

    Solution type for standard FEM.

"""
struct Solution_FEM <: AbstractSolution
    # Solution at nodes
    u :: Array{Float64,2}
    uGrad :: Array{Array{Float64,1},2}

end # end type


# --------------------------------------------------------------------
# --------------------------------------------------------------------


# Outer constructor
function Solution_FEM(u :: Array{Float64,2}, m :: Mesh.TriangleMesh.TriMesh)

    uGrad = [[zeros(2) for i in 1:m.n_cell] for j in 1:size(u,2)]
    uGrad = hcat(uGrad...)

    return Solution_FEM(u, uGrad)
end


function Solution_FEM(dof :: AbstractDof,
                          par :: Parameter.Parameter_FEM)

    # Reserve memory for the solution
    u = Array{Float64,2}(undef, dof.n_node, par.n_steps+1)
    uGrad = [[zeros(2) for i in 1:dof.n_elem] for j in 1:(par.n_steps+1)]
    uGrad = hcat(uGrad...)
    
    return Solution_FEM(u, uGrad)
end


"""
    Solution_FEM(u_in :: Array{Array{Array{Float64,1},1},1},
                        n_node :: Int64)

    Outer constructor for struct 'Solution_FEM' from post-processed FEM data
    (solution evaluated at a suitable set of points).

"""
function Solution_FEM(u_in :: Array{Array{Array{Float64,1},1},1},
                        m :: Mesh.TriangleMesh.TriMesh)

    # Reserve memory for the solution
    u = Array{Float64,2}(undef, m.n_point, length(u_in[1][1][:]))

    for i in 1:length(u_in)
        u[i,:] = u_in[i][1][:]
    end

    return Solution_FEM(u, m)
end


"""
    Solution_FEM(u_in :: Array{Array{Array{Array{Float64,1},1},1},1},
                        n_node :: Int64)

    Outer constructor for struct 'Solution_FEM' from post-processed MsFEM data
    (solution evaluated at a suitable set of points).

"""
function Solution_FEM(u_in :: Array{Array{Array{Array{Float64,1},1},1},1},
                        m :: Mesh.TriangleMesh.TriMesh)

    # Reserve memory for the solution
    u = Array{Float64,2}(undef, m.n_point, length(u_in[1][1][1][:]))

    for i in 1:length(u_in)
        u[i,:] = u_in[i][1][1][:]
    end

    return Solution_FEM(u, m)
end