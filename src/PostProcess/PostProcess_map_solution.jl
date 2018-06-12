"""
	map_solution(sol :: FEM.Solution_FEM,
						mesh :: TriangleMesh.TriMesh,
						mesh_target :: TriangleMesh.TriMesh)

	Map FEM solution to target_mesh.

"""
function map_solution(sol :: FEM.Solution_FEM,
						mesh :: TriangleMesh.TriMesh,
						mesh_target :: TriangleMesh.TriMesh)

	# We need the solution at the new mesh points
	x = mesh_target.point

	# Create an evaluation struct
	E = Evaluate_FEM(sol, mesh);

	# Evaluate the solution at the new points
	u = PostProcess.evaluate(E, mesh, x);

	# Construct a new solution on mesh_target
	sol_target = FEM.Solution_FEM(u, mesh_target)

	PostProcess.evaluateGrad!(sol_target, mesh_target)

	return sol_target
end


"""
	map_solution(sol :: FEM.Solution_MsFEM,
						mesh_collection :: Mesh.TriMesh_collection,
						mesh_target :: TriangleMesh.TriMesh)
						
	Map MsFEM solution to target_mesh. New data structure is
	'FEM.Solution_FEM'.

"""
function map_solution(sol :: FEM.Solution_MsFEM,
						mesh_collection :: Mesh.TriMesh_collection,
						mesh_target :: TriangleMesh.TriMesh)
	
	# We need the solution at the new mesh points
	x = mesh_target.point
	
	# Create an evaluation struct
	E = Evaluate_MsFEM(sol, mesh_collection);

	# Evaluate the solution at the new points
	u = PostProcess.evaluate(E, mesh_collection, x);

	# Construct a new solution on mesh_target
	# sol_target = FEM.Solution_FEM(u, mesh_target.n_point)
	sol_target = FEM.Solution_FEM(u, mesh_target)

	PostProcess.evaluateGrad!(sol_target, mesh_target)

	return sol_target
end