# ------------------------------------------------------------------------------
# -------   Implicit Euler   -------
struct ImplEuler <: AbstractTimeIntegrator

    systemData :: ImplEulerData
    dt :: Float64
    
end
# ------------------------------------------------------------------------------


function ImplEuler(dof :: FEM.AbstractDof, 
                    mesh :: Mesh.TriangleMesh.TriMesh,
                    problem :: Problem.AbstractProblem)

    systemData = ImplEulerData(dof, mesh, problem)
        
    return new(systemData)
end


# -------   System data   -------
include("TimeIntegrator_implEuler_systemData_ADE.jl")


# -------   Step function   -------
include("TimeIntegrator_implEuler_systemStep.jl")


# # -------   System update   -------
# include("TimeIntegrator_implEuler_systemUpdate_ADE.jl")
