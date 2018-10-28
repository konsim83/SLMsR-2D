struct Dof_collection{FEM_order_f} <: AbstractDofCollection

    dof :: AbstractDof
    dof_f :: Array{AbstractDof,1}

end # end type


function Dof_collection(mesh_collection :: Mesh.TriMesh_collection,
                        dof :: AbstractDof,
                        problem_f :: Array{Problem.AbstractBasisProblem,1},
		                FEM_order_f :: Int)

    dof_f = Array{Dof_Pk{FEM_order_f},1}(undef, length(mesh_collection.mesh_f))
    for i in 1:length(mesh_collection.mesh_f)
        dof_f[i] = Dof_Pk(mesh_collection.mesh_f[i], problem_f[i], FEM_order_f)
    end

    return Dof_collection{FEM_order_f}(dof, dof_f)
end # end function