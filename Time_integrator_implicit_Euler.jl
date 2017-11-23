type Implicit_Euler <: Time_stepper

    info :: String
    
    dt :: Float64

    make_step! :: Function
    
    function Implicit_Euler(par :: Parameter)
        this = new()

        this.info = "Implicit Euler method."
        
        this.dt = par.dt

        # ----------------------------------------------------------------------
        this.make_step! = function(dof :: FEM.Dof,
                                  system_data :: FEM.Setup_system_ADE_implEuler,
                                  solution :: FEM.Solution)

            """

            Accept as input an appropriately set up system in form of
            a type and do time step from k_time -> k_time+1.

            """

            solution.u[dof.ind_node_non_dirichlet,system_data.k_time+1] = dof.map_ind_dof2mesh(system_data.system_matrix \ system_data.rhs)
            solution.u[dof.ind_node_dirichlet,system_data.k_time+1] = solution.u[dof.ind_node_non_dirichlet,system_data.k_time]
            
            return nothing
        end
        # ----------------------------------------------------------------------
        
        return this
    end
    
end
