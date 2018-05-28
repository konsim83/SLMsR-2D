# ----------------------------------
# Step for a physical problem.
function makeStep!(timeInt :: ImplEuler, 
					dof :: FEM.AbstractDof,
					uNew :: T, 
					uOld :: Array{Float64}) where {T <: ArrayViewUnion}

    innd = timeInt.systemData.ind_node_non_dirichlet
    ind = timeInt.systemData.ind_node_dirichlet

    uOldDof = FEM.map_vec_mesh2dof(dof, uOld)
    uNewDof = FEM.map_vec_mesh2dof(dof, uOld)

    # Solve the system and map to mesh variables
    uNewDof[innd,:] = timeInt.systemData.system_matrix \ timeInt.systemData.system_rhs
    uNewDof[ind,:] = uOldDof[ind,:]
    
    uNew[:,:] = FEM.map_vec_dof2mesh(dof, uNewDof)

    return nothing
end
# ----------------------------------