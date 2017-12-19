type Solution_FEM <: AbstractSolution
    """
    Solution type for standard FEM.
    """

    # Solution at nodes
    u :: Array{Float64,2}

    function Solution_FEM(dof :: AbstractDof,
                          par :: Parameter_FEM)

        # Reserve memory for the solution
        u = Array{Float64,2}(dof.n_node, par.n_steps+1)
        
        return new(u)
    end
end # end type


# ----------------------------------------------------------------------------------------


type Solution_MsFEM <: AbstractSolution
    """
    Solution type for multiscale FEM.
    """

    # Solution at nodes
    u :: Array{Float64,2}

    # Multiscale basis functions
    phi_1 :: Array{Array{Float64,2},1}
    phi_2 :: Array{Array{Float64,2},1}
    phi_3 :: Array{Array{Float64,2},1}

    # Storage of element matrices
    mass :: Array{SparseMatrixCSC{Float64,Int64},2}
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
        
        mass = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        advection = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        diffusion = Array{SparseMatrixCSC{Float64,Int64},2}(dof.dof.n_elem, par.n_steps+1)
        
        # Reserve memory for each element of the uninitialized array
        for i=1:dof.dof.n_elem
            phi_1[i] = Array{Float64,2}(dof.dof_f[i].n_node, par.n_steps+1)
            phi_2[i] = Array{Float64,2}(dof.dof_f[i].n_node, par.n_steps+1)
            phi_3[i] = Array{Float64,2}(dof.dof_f[i].n_node, par.n_steps+1)

            pattern_matrix = sparse(dof.dof_f[i].ind_test, dof.dof_f[i].ind, zeros(Float64, length(dof.dof_f[i].ind)), dof.dof_f[i].n_true_dof, dof.dof_f[i].n_true_dof)
            for j=1:par.n_steps+1
                mass[i,j] = copy(pattern_matrix)
                advection[i,j] = copy(pattern_matrix)
                diffusion[i,j] = copy(pattern_matrix)
            end
        end
        
        return new(u, phi_1, phi_2, phi_3, mass, advection, diffusion)
    end
end # end type
