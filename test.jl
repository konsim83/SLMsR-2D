if true
	include("Parameter.jl")
	include("Geometry.jl")
	include("Problem.jl")
	include("Quad.jl")
	include("Mesh.jl")
	include("FiniteDiff.jl")
	include("FEM.jl")
    include("Time_integrator.jl")
    include("Solver.jl")
    include("PostProcess.jl")
    include("Vis.jl")
end


periodic = false

if periodic
	par_global = Parameter.Parameter_FEM(1.0, 1/500, 2, 0, 1, 3, 1)
	m = Mesh.mesh_unit_square(5)
	r = FEM.RefEl_Pk(1)
	p = Problem.Gaussian(1.0)
	tri = Geometry.Triangle(r.node)
	pb = Problem.BasisFun(p, tri)
	per = Mesh.DoublePeriodicUnitSquare()
	q = Quad.Quad_simplex(6)
	d = FEM.Dof_Pk_periodic(m, p, per, 1)
	ts = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(d, m, p)
else
	par_global = Parameter.Parameter_FEM(1.0, 1/500, 2, 0, 1, 3, 1)
	m = Mesh.mesh_unit_simplex(5)
	r = FEM.RefEl_Pk(1)
	p = Problem.Gaussian(1.0)
	tri = Geometry.Triangle(r.node)
	pb = Problem.BasisFun(p, tri)

	q = Quad.Quad_simplex(6)
	d = FEM.Dof_Pk(m, pb, 1)
	ts = Time_integrator.ImplEuler{Time_integrator.System_data_implEuler_ADE}(d, m, p)
end