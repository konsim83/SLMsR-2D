module TimeIntegrator

using FEM, Quad, Parameter, Mesh, Problem


abstract type AbstractTimeIntegrator end
abstract type AbstractSystemData end


# ------------------------------------------------------------------------------
# -------   ImplEuler   -------
include("TimeIntegrator_implEuler.jl")

end # end module