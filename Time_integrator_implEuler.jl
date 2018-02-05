# ------------------------------------------------------------------------------
# -------   Implicit Euler   -------
type ImplEuler{T <: AbstractSystem_data_implEuler} <: AbstractTime_integrator

    system_data :: T
    
    function ImplEuler{T}(dof :: FEM.AbstractDof, mesh :: Mesh.TriangleMesh.TriMesh, problem :: Problem.AbstractProblem) where {T <: AbstractSystem_data_implEuler}
        
	# Reserve Memory for System data only
        if problem.type_info=="ADE"
	    system_data = System_data_implEuler_ADE(dof, mesh, problem)
        else
            error("Problem type not implemented yet.")
        end
            
        return new(system_data)
    end
end
# ------------------------------------------------------------------------------




# -------   System data   -------
include("Time_integrator_implEuler_systemData_ADE.jl")


# -------   Step function   -------
include("Time_integrator_implEuler_systemStep_ADE.jl")


# -------   System update   -------
include("Time_integrator_implEuler_systemUpdate_ADE.jl")
