module TimeIntegrator

using FEM, Quad, Parameter, Mesh, Problem



abstract type AbstractTimeIntegrator end

abstract type AbstractSystem_data end
abstract type AbstractSystem_data_implEuler <: AbstractSystem_data end



# ------------------------------------------------------------------------------
# -------   ImplEuler   -------
include("TimeIntegrator_implEuler.jl")


end # end module
