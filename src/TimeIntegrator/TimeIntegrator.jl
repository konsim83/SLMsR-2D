module TimeIntegrator

using FEM, Quad, Parameter, Mesh, Problem

ArrayViewUnion = Union{ SubArray{Float64,1,Array{Float64,2},Tuple{Base.Slice{Base.OneTo{Int64}},Int64},true},
						SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}
						}

abstract type AbstractTimeIntegrator end
abstract type AbstractSystemData end


# ------------------------------------------------------------------------------
# -------   ImplEuler   -------
include("TimeIntegrator_implEuler.jl")

end # end module