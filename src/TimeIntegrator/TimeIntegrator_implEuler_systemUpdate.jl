function updateSystem!(systemData :: ImplEulerData,
                        M :: SparseMatrixCSC{Float64,Int64},
                        A :: SparseMatrixCSC{Float64,Int64},
                        f :: Array{Float64},
                        dof :: FEM.AbstractDof,
                        uOld :: Array{Float64},
                        dt :: Float64)
    
    # Just an abbreviation
    innd = systemData.ind_node_non_dirichlet
    ind = systemData.ind_node_dirichlet

    uOldDof = FEM.map_vec_mesh2dof(dof, uOld[:,:])
    fDof = FEM.map_vec_mesh2dof(dof, f[:,:])

    systemData.system_matrix[:,:] = (M - dt*A)[innd,innd]
    
    systemData.system_rhs[:,:] = M[innd,innd]*uOldDof[innd,:] + dt*fDof[innd,:]
    
    if !isempty(ind)
        systemData.system_rhs[:,:] +=  M[innd,ind]*uOldDof[ind,:] - (M - dt*A)[innd,ind]*uOldDof[ind,:]
    end
    
    return nothing
end


function updateSystem!(systemData :: ImplEulerData,
                        M :: SparseMatrixCSC{Float64,Int64},
                        A :: SparseMatrixCSC{Float64,Int64},
                        f :: Array{Float64},
                        dof :: FEM.AbstractDof,
                        uOld :: Array{Float64},
                        uNewBC :: Array{Float64},
                        dt :: Float64)
    
    # Just an abbreviation
    innd = systemData.ind_node_non_dirichlet
    ind = systemData.ind_node_dirichlet

    uOldDof = FEM.map_vec_mesh2dof(dof, uOld[:,:])
    uNewBCDof = FEM.map_vec_mesh2dof(dof, uNewBC[:,:])
    fDof = FEM.map_vec_mesh2dof(dof, f[:,:])

    systemData.system_matrix[:,:] = (M - dt*A)[innd,innd]
    
    systemData.system_rhs[:,:] = M[innd,innd]*uOldDof[innd,:] + dt*fDof[innd,:]
    
    if !isempty(ind)
        systemData.system_rhs[:,:] +=  M[innd,ind]*uOldDof[ind,:] - (M - dt*A)[innd,ind]*uNewBCDof[ind,:]
    end
    
    return nothing
end