function SemiLagrange(solution :: FEM.Solution_MsFEM,
						mesh_collection :: Mesh.TriMesh_collection,
                        dof_collection :: FEM.Dof_collection,
                        par :: Parameter.Parameter_MsFEM,
                        problem :: Problem.AbstractPhysicalProblem,
                        problem_f :: Array{Problem.AbstractBasisProblem,1})

	# Needs to be set up only once since velocity is the same everywhere
	velocity = function(x :: Array{Float64,2}, parameter, t :: Float64)

		return hcat(-Problem.velocity(problem, -t, x)...)
	end

	for i in 1:mesh_collection.mesh.n_cell
		point = mesh_collection.mesh_f[i].point
		point_orig = traceback(point, T, velocity, par)
	end


end


function traceback(point :: Array{Float64,2}, T :: Float64, 
					velocity :: Function, par :: Parameter.Parameter_MsFEM,)

	tspan = (-T, -T+par.dt)
	prob = ODEProblem(velocity, point, tspan)
	sol = solve(prob)

	return point_orig
end