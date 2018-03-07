# ----------------------------------
# Step for a physical problem.
function make_step!(solution :: FEM.AbstractSolution,
                    time_int :: ImplEuler,
                    mesh :: Mesh.TriangleMesh.TriMesh,
                    dof :: FEM.AbstractDof,
                    ref_el :: FEM.AbstractRefEl,
                    quad :: Quad.AbstractQuad,
                    par :: Parameter.AbstractParameter,
                    problem :: Problem.AbstractPhysicalProblem,
                    k_time :: Int64)
    

    """

    Accept as input an Implicit Euler system in form of a type
    and do time step from k_time -> k_time+1.

    """


    # ---------------------------------------
    # Manipulate the system data locally
    update_system!(solution,
                   time_int.system_data,
                   mesh,
                   dof,
                   ref_el,
                   quad,
                   par,
                   problem,
                   k_time)
    # ---------------------------------------

    # Solve the system and map to mesh variables
    time_int.system_data.u_temp[time_int.system_data.ind_node_non_dirichlet] = time_int.system_data.system_matrix \ time_int.system_data.system_rhs
    solution.u[:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp)
    
    return nothing
end
# ----------------------------------


# ----------------------------------
# Step for a basis problem.
function make_step!(solution :: FEM.AbstractSolution,
                    time_int :: ImplEuler,
                    mesh :: Mesh.TriangleMesh.TriMesh,
                    dof :: FEM.AbstractDof,
                    ref_el :: FEM.AbstractRefEl,
                    quad :: Quad.AbstractQuad,
                    par :: Parameter.AbstractParameter,
                    problem :: Problem.AbstractBasisProblem,
                    k_time :: Int64,
                    ind_cell :: Int64)
    

    """

    Accept as input an Implicit Euler system in form of a type
    and do time step from k_time -> k_time+1.

    """


    # ---------------------------------------
    # Manipulate the system data locally
    update_system!(solution,
                   time_int.system_data,
                   mesh,
                   dof,
                   ref_el,
                   quad,
                   par,
                   problem,
                   k_time,
                   ind_cell)
    # ---------------------------------------

    
    # Solve the system and map to mesh variables
    time_int.system_data.u_temp[time_int.system_data.ind_node_non_dirichlet,:] = time_int.system_data.system_matrix \ time_int.system_data.system_rhs
    solution.phi_1[ind_cell][:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp[:,1])
    solution.phi_2[ind_cell][:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp[:,2])
    solution.phi_3[ind_cell][:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp[:,3])    
    
    return nothing
end
# ----------------------------------
