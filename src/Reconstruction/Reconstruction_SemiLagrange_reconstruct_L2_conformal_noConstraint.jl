function reconstruct_L2_noConstraint_conformal(solution :: FEM.Solution_MsFEM,
                                                mesh_collection :: Mesh.TriMesh_collection,
                                                par :: Parameter.Parameter_MsFEM,
                                                problem :: Problem.AbstractPhysicalProblem,
                                                problem_f :: Problem.AbstractBasisProblem,
                                                point :: Array{Float64,2},
                                                ind_cell :: Int)
    
    # ------------------------------------------------
    # Evaluate the solution at the traced back points
    u_orig = Problem.u_init(problem, point)

    uOpt = reconstruct_L2_noConstraint_conformal(solution,
                                                    mesh_collection,
                                                    par,
                                                    problem_f,
                                                    u_orig,
                                                    2,
                                                    ind_cell)
    
    return uOpt
end


function reconstruct_L2_noConstraint_conformal(solution :: FEM.Solution_MsFEM,
                                                mesh_collection :: Mesh.TriMesh_collection,
                                                par :: Parameter.Parameter_MsFEM,
                                                problem_f :: Problem.AbstractBasisProblem,
                                                u_orig :: Array{Float64,1},
                                                k_time :: Int,
                                                ind_cell :: Int)

    # ------------------------------------------------
    # Evaluate the solution at the traced back points
    # u = PostProcess.evaluate(solution, mesh_collection, point_orig, k_time-1)
    u = u_orig
    
	ind_corner = [sort(m_f.segment[:,m_f.segment_marker.==1],dims=1)[1] ; 
                    sort(m_f.segment[:,m_f.segment_marker.==1],dims=1)[end] ; 
                    sort(m_f.segment[:,m_f.segment_marker.==2],dims=1)[end] ]
    # U = solution.u[mesh_collection.mesh.cell[:,ind_cell],k_time-1]
    U = u_orig[ind_corner]

    u1 = U[1]
    u2 = U[2]
    u3 = U[3]
    # ------------------------------------------------

    u0 = Problem.u_init(problem_f, mesh_collection.mesh_f[ind_cell].point)

    k = par.k
    m_f = mesh_collection.mesh_f[ind_cell]
    n_dof = m_f.n_point
    uOpt = zeros(n_dof,3)


    # ------------------------------------------------
    # Reconstruct edges
    reconstruct_edge_L2_noConstraint!(uOpt, m_f, U, u0, u, 1, par.k[4])
    reconstruct_edge_L2_noConstraint!(uOpt, m_f, U, u0, u, 2, par.k[4])
    reconstruct_edge_L2_noConstraint!(uOpt, m_f, U, u0, u, 3, par.k[4])
    # ------------------------------------------------

    # ------------------------------------------------
    # Reconstruct interior
    uOpt = vec(uOpt)

    # Contraints for nodal values of basis at boundary for ALL THREE BASIS
    # FUNCTIONS!!!
    ind_con = [findall(vec(m_f.point_marker).!=0) ; 
                findall(vec(m_f.point_marker).!=0).+n_dof ;
                findall(vec(m_f.point_marker).!=0).+2*n_dof]

    constr_val = uOpt[ind_con]
    ind_uncon = setdiff(1:(3*n_dof), ind_con)

    n = length(u)
    system_matrix = [(u1^2+k[1])*speye(n) u1*u2*speye(n) u1*u3*speye(n) ; 
                        u2*u1*speye(n) (u2^2+k[2])*speye(n) u2*u3*speye(n) ; 
                        u3*u1*speye(n) u3*u2*speye(n) (u3^2+k[3])*speye(n)]
    rhs = [k[1]*u0[:,1] + u1*u; 
            k[2]*u0[:,2] + u2*u ; 
            k[3]*u0[:,3] + u3*u]#  - system_matrix[:,ind_con]*constr_val

    # uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]
    uOpt = system_matrix \ rhs
    # ------------------------------------------------

    return reshape(uOpt,:,3)
end
