type Implicit_Euler <: Time_stepper

    info :: String
    
    dt :: Float64

    make_step :: Function
    
    function Implicit_Euler(par :: Parameter)
        this = new()

        this.info = "Implicit Euler method."
        
        this.dt = par.dt

        # ----------------------------------------------------------------------
        this.make_step = function(solution :: FEM.Solution,
                                  mesh :: Mesh.TriMesh,
                                  dof :: FEM.Dof,
                                  ref_el :: FEM.RefEl,
                                  par :: Parameter.Parameter_top,
                                  problem :: Problem.Problem_top,
                                  k_time :: Int64)

            """
            Make a time step from k_time to k_time+1 by mutation
            of fields in the solution variable.
            """

            M = FEM.assemble_mass(solution,
                                  mesh,
                                  dof,
                                  ref_el,
                                  par,
                                  problem,
                                  k_time+1)
            
            A = FEM.assemble_advection(solution,
                                       mesh,
                                       dof,
                                       ref_el,
                                       par,
                                       problem,
                                       k_time+1)

            D = FEM.assemble_diffusion(solution,
                                       mesh,
                                       dof,
                                       ref_el,
                                       par,
                                       problem,
                                       k_time+1)

            if dof.n_node_neumann > 0
                error("Neumann boundary integral not implemented yet.")
            end

            ind_node_non_dirichlet[]
            
            #solve_system(this.dt, M, A, D, solution, )
        end
        # ----------------------------------------------------------------------
        
        return this
    end
    
end
