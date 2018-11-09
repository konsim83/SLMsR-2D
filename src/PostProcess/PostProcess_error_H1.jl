function error_H1(s1 :: FEM.Solution_FEM,
					s2 :: FEM.Solution_FEM)

	size(s1.uGrad)!=size(s2.uGrad) ? error("Solutions must be of same size.") :

	diffGrad = s2.uGrad - s1.uGrad

	err = sqrt.(sum([sum(a.^2) for a in diffGrad],dims=1) ./ sum([sum(a.^2) for a in s1.uGrad],dims=1))

	return err
end