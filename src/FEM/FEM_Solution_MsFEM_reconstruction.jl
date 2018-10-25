"""
    struct Solution_MsFEM_reconstruction <: AbstractSolution

    Solution type for multiscale FEM.

"""
struct Solution_MsFEM_reconstruction{FEM_order} <: AbstractSolution
    # Solution at nodes
    u :: Array{Float64,2}

    # Standard basis functions
    phi :: Array{Array{Array{Float64,2},1},1}

    # Time derivatives of standard basis functions
    phi_t :: Array{Array{Array{Float64,2},1},1}
    
    # Storage of element matrices
    # mass :: Array{SparseMatrixCSC{Float64,Int64},1}
    # advection :: Array{SparseMatrixCSC{Float64,Int64},2}
    # diffusion :: Array{SparseMatrixCSC{Float64,Int64},2}
    
end # end type



function Solution_MsFEM_reconstruction(dof :: Dof_collection,
                                            par :: Parameter_MsFEM)
    if FEM_order==1
        # --------------------------------------------
        # Reserve memory for the solution
        u = Array{Float64, 2}(undef, dof.dof.n_node, par.n_steps+1)
        
        # Set up an array of arrays for the basis
        phi = [[Array{Float64,2}(undef, dof.dof_f[i].n_true_dof, par.n_steps+1) for i in 1:dof.dof.n_elem] for j in 1:3]


        # Set up an array of arrays for the basis
        phi_t = [[Array{Float64,2}(undef, dof.dof_f[i].n_true_dof, par.n_steps+1) for i in 1:dof.dof.n_elem] for j in 1:3]
        
        # mass = Array{SparseMatrixCSC{Float64,Int64},1}(dof.dof.n_elem)
        # advection = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        # diffusion = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        
        # # Reserve memory for each element of the uninitialized array
        # for i=1:dof.dof.n_elem
        #     pattern_matrix = sparse(dof.dof_f[i].ind_test, dof.dof_f[i].ind, 
        #                             zeros(Float64, length(dof.dof_f[i].ind)), 
        #                             dof.dof_f[i].n_true_dof, dof.dof_f[i].n_true_dof)
        #     mass[i] = copy(pattern_matrix)
        #     for j=1:par.n_steps+1
        #         advection[i,j] = copy(pattern_matrix)
        #         diffusion[i,j] = copy(pattern_matrix)
        #     end
        # end
        # --------------------------------------------
    else
        # --------------------------------------------
        error("FEM order not implemented for reconstruction.")
        # --------------------------------------------
    end
    return Solution_MsFEM_reconstruction{FEM_order}(u, phi, phi_t)
end