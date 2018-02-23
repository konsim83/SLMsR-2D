"""
	error_L2(sol1 :: FEM.Solution_FEM, sol2 :: FEM.Solution_FEM)

	Compute the L2-difference of two FEM solutions on the same grid.

"""
function error_L2(sol1 :: FEM.Solution_FEM, sol2 :: FEM.Solution_FEM)

	size(sol1.u)!=size(sol2.u) ? error("Solutions must be of same size.") :

	err = sqrt.(sum(abs.(sol1.u-sol2.u).^2,1) ./ sum(abs.(sol1.u).^2,1))

	return err
end


"""
	error_max(sol1 :: FEM.Solution_FEM, sol2 :: FEM.Solution_FEM)

	Compute the maximum norm of the difference of two FEM solutions on the
	same grid.

"""
function error_max(sol1 :: FEM.Solution_FEM, sol2 :: FEM.Solution_FEM)

	size(sol1.u)!=size(sol2.u) ? error("Solutions must be of same size.") :

	err = (maximum(abs.(sol1.u-sol2.u).^2,1) ./ maximum(abs.(sol1.u).^2,1))

	return err
end