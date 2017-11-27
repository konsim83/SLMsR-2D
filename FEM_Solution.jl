type Solution_FEM <: Solution
    """
    Solution type for standard FEM.
    """

    # Solution at nodes
    u :: Array{Float64,2}

    function Solution_FEM(dof :: Dof,
                          par :: Parameter_FEM)
        
        this = new()

        # Reserve memory for the solution
        this.u = Array(Float64, dof.n_node, par.n_steps)
        
        return this
    end
end # end type


# ----------------------------------------------------------------------------------------


type Solution_MsFEM <: Solution
    """
    Solution type for multiscale FEM.
    """

    # Solution at nodes
    u :: Array{Float64,2}

    # Multiscale basis functions
    phi_1 :: Array{Array{Float64,2},1}
    phi_2 :: Array{Array{Float64,2},1}
    phi_3 :: Array{Array{Float64,2},1}

    function Solution_MsFEM(dof :: Dof_square,
                          par :: Parameter_MsFEM)
        
        this = new()

        # Reserve memory for the solution
        this.u = Array(Float64, dof.n_node, par.n_steps)

        # Set up an array of arrays for the basis
        this.phi_1 = Array(Array{Float64,2}, par.n_elem)
        this.phi_2 = Array(Array{Float64,2}, par.n_elem)
        this.phi_3 = Array(Array{Float64,2}, par.n_elem)

        # Reserve memory for each element of the uninitialized array
        for i=1:par.n_elem
            this.phi_1[i] = Array(Float64, dof.n_nodes, par.n_steps)
            this.phi_2[i] = Array(Float64, dof.n_nodes, par.n_steps)
            this.phi_3[i] = Array(Float64, dof.n_nodes, par.n_steps)
        end
        
        return this
    end
end # end type
