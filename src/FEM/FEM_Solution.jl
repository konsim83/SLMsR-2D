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


# Outer constructor
function Solution_FEM(dof :: AbstractDof,
                          par :: Parameter_FEM)

    # Reserve memory for the solution
    u = Array{Float64,2}(dof.n_node, par.n_steps+1)
    
    return Solution_FEM(u)
end


# ----------------------------------------------------------------------------------------


"""
    struct Solution_MsFEM <: AbstractSolution

    Solution type for multiscale FEM.

"""
struct Solution_MsFEM <: AbstractSolution
    # Solution at nodes
    u :: Array{Float64,2}

    # Multiscale basis functions
    phi_1 :: Array{Array{Float64,2},1}
    phi_2 :: Array{Array{Float64,2},1}
    phi_3 :: Array{Array{Float64,2},1}

    # Time derivatives of multiscale basis functions
    phi_1_t :: Array{Array{Float64,2},1}
    phi_2_t :: Array{Array{Float64,2},1}
    phi_3_t :: Array{Array{Float64,2},1}

    # Storage of element matrices
    mass :: Array{SparseMatrixCSC{Float64,Int64},1}
    advection :: Array{SparseMatrixCSC{Float64,Int64},2}
    diffusion :: Array{SparseMatrixCSC{Float64,Int64},2}
    
    function Solution_MsFEM(dof :: Dof_collection,
                            par :: Parameter_MsFEM)
        
        # Reserve memory for the solution
        u = Array{Float64, 2}(dof.dof.n_node, par.n_steps+1)
        
        # Set up an array of arrays for the basis
        phi_1 = Array{Array{Float64,2}}(dof.dof.n_elem)
        phi_2 = Array{Array{Float64,2}}(dof.dof.n_elem)
        phi_3 = Array{Array{Float64,2}}(dof.dof.n_elem)

        # Set up an array of arrays for the basis
        phi_1_t = Array{Array{Float64,2}}(dof.dof.n_elem)
        phi_2_t = Array{Array{Float64,2}}(dof.dof.n_elem)
        phi_3_t = Array{Array{Float64,2}}(dof.dof.n_elem)
        
        mass = Array{SparseMatrixCSC{Float64,Int64},1}(dof.dof.n_elem)
        advection = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        diffusion = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        
        # Reserve memory for each element of the uninitialized array
        for i=1:dof.dof.n_elem
            phi_1[i] = Array{Float64,2}(dof.dof_f[i].n_true_dof, par.n_steps+1)
            phi_2[i] = Array{Float64,2}(dof.dof_f[i].n_true_dof, par.n_steps+1)
            phi_3[i] = Array{Float64,2}(dof.dof_f[i].n_true_dof, par.n_steps+1)

            phi_1_t[i] = Array{Float64,2}(dof.dof_f[i].n_true_dof, par.n_steps+1)
            phi_2_t[i] = Array{Float64,2}(dof.dof_f[i].n_true_dof, par.n_steps+1)
            phi_3_t[i] = Array{Float64,2}(dof.dof_f[i].n_true_dof, par.n_steps+1)

            pattern_matrix = sparse(dof.dof_f[i].ind_test, dof.dof_f[i].ind, zeros(Float64, length(dof.dof_f[i].ind)), dof.dof_f[i].n_true_dof, dof.dof_f[i].n_true_dof)
            mass[i] = copy(pattern_matrix)
            for j=1:par.n_steps+1
                advection[i,j] = copy(pattern_matrix)
                diffusion[i,j] = copy(pattern_matrix)
            end
        end
        
        return new(u, phi_1, phi_2, phi_3, phi_1_t, phi_2_t, phi_3_t, mass, advection, diffusion)
    end
end # end type
