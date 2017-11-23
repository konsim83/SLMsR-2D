"""

This is the part of the system that you need to program yourself. Its
purpose is to setup the system correctly that can then be passed to a
generic time integrator. The data that is produced here needs to fit
the integrator.

"""
type Setup_system_ADE_implEuler <: System

    system_matrix :: Array{Float64,2}
    rhs :: Array{Float64,1}

    k_time :: Int64

    function Setup_system_ADE_implEuler(solution :: FEM.Solution,
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
        
        this = new()

        this.k_time = k_time
        
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

        this.system_matrix = (M-par.dt*(D-A))[dof.ind_nodes_non_dirichlet,dof.ind_nodes_non_dirichlet]      

        if dof.is_periodic
            this.rhs = M[dof.ind_nodes_non_dirichlet,dof.ind_nodes_non_dirichlet] * dof.map_ind_mesh2_dof(solution.u[:,k_time])
        else
            
            this.rhs = (M[dof.ind_nodes_non_dirichlet,dof.ind_nodes_non_dirichlet] * solution.u[dof.ind_nodes_non_dirichlet,k_time] -
                        M[dof.ind_nodes_non_dirichlet,dof.ind_nodes_dirichlet] * solution.u[dof.ind_nodes_dirichlet,k_time])
        end
        
        return this
    end
    
end
