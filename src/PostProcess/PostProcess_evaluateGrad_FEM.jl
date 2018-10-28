
function evaluateGrad!(s :: FEM.Solution_FEM, m :: TriangleMesh.TriMesh)

	uGrad = [[zeros(2) for i in 1:m.n_cell] for j in 1:size(s.u,2)]
	uGrad = hcat(uGrad...)

	for i in 1:m.n_cell
		ind_p = m.cell[:,i]
		P_inv = transpose([m.point[:,ind_p] ; ones(1,3)]) \ I
		for j in 1:size(s.u,2)
			basis_weight = P_inv * diagm(s.u[ind_p,j])
			uGrad[i,j] = sum(basis_weight[1:2,:], 2)[:]
		end
	end

	s.uGrad[:,:] = uGrad[:,:]
end