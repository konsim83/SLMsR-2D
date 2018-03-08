# --------------------------------------------------
# --------------------------------------------------

# path = pwd()
# println("Adding to LOAD_PATH:   $path")
# push!(LOAD_PATH, pwd())

# for (root, dirs, files) in walkdir(pwd()*"/src")
#     for dir in dirs
#     	path = joinpath(root, dir)
#     	println("Adding to LOAD_PATH:   $path")
#         push!(LOAD_PATH, path)
#     end
# end

# --------------------------------------------------
# --------------------------------------------------


using Mesh, Parameter, Problem_Test, Quad, FiniteDiff, Geometry, FEM

using Base.Test

# include("Test_Polygon.jl")

include("Test_assembly_std.jl")