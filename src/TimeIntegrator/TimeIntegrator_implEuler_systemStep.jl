# ----------------------------------
# Step for a physical problem.
function makeStep!(timeInt :: ImplEuler, uNew :: Array{Float64}, uOld :: Array{Float64})

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

    uOldDof = FEM.map_vec_mesh2dof(dof, uOld)
    uNewDof = FEM.map_vec_mesh2dof(dof, uOld)

    # Solve the system and map to mesh variables
    uNewDof[timeInt.systemData.ind_node_non_dirichlet] = timeInt.system_data.system_matrix \ timeInt.system_data.system_rhs
    uNew[:,:] = FEM.map_vec_dof2mesh(dof, uNewDof)
    
    return nothing
end
# ----------------------------------

# # ----------------------------------
# # Step for a basis problem.
# function make_step!(solution :: FEM.AbstractSolution,
#                     time_int :: ImplEuler,
#                     mesh :: Mesh.TriangleMesh.TriMesh,
#                     dof :: FEM.AbstractDof,
#                     ref_el :: FEM.AbstractRefEl,
#                     quad :: Quad.AbstractQuad,
#                     par :: Parameter.AbstractParameter,
#                     problem :: Problem.AbstractBasisProblem,
#                     k_time :: Int64,
#                     ind_cell :: Int64)

#     # ---------------------------------------
#     # Manipulate the system data locally
#     update_system!(solution,
#                    time_int.system_data,
#                    mesh,
#                    dof,
#                    ref_el,
#                    quad,
#                    par,
#                    problem,
#                    k_time,
#                    ind_cell)
#     # ---------------------------------------

    
#     # Solve the system and map to mesh variables
#     time_int.system_data.u_temp[time_int.system_data.ind_node_non_dirichlet,:] = time_int.system_data.system_matrix \ time_int.system_data.system_rhs
#     solution.phi_1[ind_cell][:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp[:,1])
#     solution.phi_2[ind_cell][:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp[:,2])
#     solution.phi_3[ind_cell][:,k_time+1] = FEM.map_vec_dof2mesh(dof, time_int.system_data.u_temp[:,3])    
    
#     return nothing
# end
# # ----------------------------------