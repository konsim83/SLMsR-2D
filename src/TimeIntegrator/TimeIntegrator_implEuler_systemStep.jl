# ----------------------------------
# Step for a physical problem.
function makeStep!(timeInt :: ImplEuler, uNew :: Array{Float64}, uOld :: Array{Float64})

    innd = timeInt.systemData.ind_node_non_dirichlet

    uOldDof = FEM.map_vec_mesh2dof(dof, uOld)
    uNewDof = FEM.map_vec_mesh2dof(dof, uOld)

    # Solve the system and map to mesh variables
    uNewDof[innd,:] = timeInt.system_data.system_matrix \ timeInt.system_data.system_rhs
    uNew[:,:] = FEM.map_vec_dof2mesh(dof, uNewDof)
    
    return nothing
end
# ----------------------------------