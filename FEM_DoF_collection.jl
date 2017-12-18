type Dof_collection{FEM_order_f}

    dof :: AbstractDof
    dof_f :: Array{AbstractDof,1}

    function Dof_collection{FEM_order_f}(mesh_collection :: Mesh.TriMesh_collection) where {FEM_order_f}

        dof = Dof_Pk_periodic_square{1}(mesh_collection.mesh)

        dof_f = Array{Dof_Pk{FEM_order_f},1}(length(mesh_collection.mesh_f))
        for i in 1:length(mesh_collection.mesh_f)
            dof_f[i] = Dof_Pk{FEM_order_f}(mesh_collection.mesh_f[i])
        end

        return new(dof, dof_f)
    end # end function    
end # end type
