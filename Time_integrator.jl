module Time_integrator

using FEM, Parameter, Mesh, Problem


# -------   Implicit Euler   -------
function make_step!(dof :: FEM.Dof,
                    system_data :: FEM.Setup_system_ADE_implEuler,
                    solution :: FEM.Solution)

            """

            Accept as input an Implicit Euler system in form of a type
            and do time step from k_time -> k_time+1.

            """

            solution.u[dof.ind_node_non_dirichlet,system_data.k_time+1] = dof.map_ind_dof2mesh(system_data.system_matrix \ system_data.rhs)
            solution.u[dof.ind_node_dirichlet,system_data.k_time+1] = solution.u[dof.ind_node_non_dirichlet,system_data.k_time]
            
    return nothing
end
# ----------------------------------

end
