module Time_integrator

using FEM, Quad, Parameter, Mesh, Problem



abstract type AbstractTime_integrator end

abstract type AbstractSystem_data end
abstract type AbstractSystem_data_implEuler <: AbstractSystem_data end



# ------------------------------------------------------------------------------
# -------   ImplEuler   -------
include("Time_integrator_implEuler.jl")


end # end module
