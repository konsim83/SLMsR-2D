# -------   System data   -------
include("TimeIntegrator_implEuler_systemData.jl")



# ----------------------------------
# -------   Implicit Euler   -------
struct ImplEuler <: AbstractTimeIntegrator

    systemData :: ImplEulerData
    
end


function ImplEuler(dof :: FEM.AbstractDof, 
                    mesh :: Mesh.TriangleMesh.TriMesh,
                    problem :: Problem.AbstractProblem)

    systemData = ImplEulerData(dof, mesh, problem)
        
    return ImplEuler(systemData)
end
# ----------------------------------
# ----------------------------------



# -------   System data   -------
include("TimeIntegrator_implEuler_systemUpdate.jl")

# -------   Step function   -------
include("TimeIntegrator_implEuler_systemStep.jl")