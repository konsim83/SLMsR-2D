@testset "FEM - Assembly Procedures - std" begin
	    
       par_global = Parameter.Parameter_FEM(1.0, 1/500, 2, 0, 1, 3, 1)
       m = Mesh.mesh_unit_square(500)
       r = FEM.RefEl_Pk(1)
       p = Problem.Gaussian(1.0)
       per = Mesh.DoublePeriodicUnitSquare()
       q = Quad.Quad_simplex(6)
       d = FEM.Dof_Pk_periodic(m, p, per, 1)

       M = FEM.assemble_mass(m,d,r,q,p)
       A = FEM.assemble_advection(m,d,r,q,par_global,p,1)
       D = FEM.assemble_diffusion(m,d,r,q,par_global,p,1)


       @test isapprox(sum(D,1), 0)
       @test isapprox(sum(D,2), 0)
       @test isapprox(maximum(abs.(D-D')), 0)
       
       @test isapprox(sum(M),1)
       @test isapprox(maximum(abs.(M-M')), 0)
       
       @test isapprox(maximum(abs.(A+A')), 0)
       @test isapprox(sum(abs.(real(eig(full(A))[1]))))
end


@testset "FEM - Assembly Procedures - std (in place)" begin
           
       par_global = Parameter.Parameter_FEM(1.0, 1/500, 2, 0, 1, 3, 1)
       m = Mesh.mesh_unit_square(500)
       r = FEM.RefEl_Pk(1)
       p = Problem.Gaussian(1.0)
       per = Mesh.DoublePeriodicUnitSquare()
       q = Quad.Quad_simplex(6)
       d = FEM.Dof_Pk_periodic(m, p, per, 1)

       M0 = FEM.assemble_mass(m,d,r,q,p)
       A0 = FEM.assemble_advection(m,d,r,q,par_global,p,1)
       D0 = FEM.assemble_diffusion(m,d,r,q,par_global,p,1)

       M = sparse(d.ind_test, d.ind, vec(mat_local), d.n_true_dof, d.n_true_dof)
       A = sparse(d.ind_test, d.ind, vec(mat_local), d.n_true_dof, d.n_true_dof)
       D = sparse(d.ind_test, d.ind, vec(mat_local), d.n_true_dof, d.n_true_dof)

       FEM.assemble_mass!(M,m,d,r,q,p)
       FEM.assemble_advection!(A,m,d,r,q,par_global,p,1)
       FEM.assemble_diffusion!(D,m,d,r,q,par_global,p,1)

       @test isapprox(M, M0)
       @test isapprox(A, A0)
       @test isapprox(D, D0)
       
end