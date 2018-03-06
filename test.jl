if true
	include("Parameter.jl")
	include("Geometry.jl")
	include("Problem.jl")
	include("Quad.jl")
	include("Mesh.jl")
	include("FiniteDiff.jl")
	include("FEM.jl")
end


par_global = Parameter.Parameter_FEM(1.0, 1/500, 2, 0, 1, 3, 1)
m = Mesh.mesh_unit_square(3)
r = FEM.RefEl_Pk(1)
p = Problem.Gaussian(1.0)
per = Mesh.DoublePeriodicUnitSquare()
q = Quad.Quad_simplex(6)
@time d = FEM.Dof_Pk_periodic(m, p, per, 1)


