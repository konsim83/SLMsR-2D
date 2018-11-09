function reconstruct_H1_conformal(solution :: FEM.Solution_MsFEM,
    									mesh_collection :: Mesh.TriMesh_collection,
                                        dof_collection :: FEM.AbstractDofCollection,
                                        ref_el_f :: FEM.AbstractRefEl,
                                        quad_f :: Quad.AbstractQuad,
    									par :: Parameter.Parameter_MsFEM,
    									problem :: Problem.AbstractPhysicalProblem,
    									problem_f :: Problem.AbstractBasisProblem,
    									point :: Array{Float64,2},
    									ind_cell :: Int)
    
    # Evaluate the solution at the traced back points
    u_orig = Problem.u_init(problem, point)

    uOpt = reconstruct_H1_conformal(solution,
                                            mesh_collection,
                                            dof_collection,
                                            ref_el_f,
                                            quad_f,
                                            par,
                                            problem_f,
                                            u_orig,
                                            2,
                                            ind_cell)
    
	return uOpt
 
end



# ===========================================================================================



function reconstruct_H1_conformal(solution :: FEM.Solution_MsFEM,
                                        mesh_collection :: Mesh.TriMesh_collection,
                                        dof_collection :: FEM.AbstractDofCollection,
                                        ref_el_f :: FEM.AbstractRefEl,
                                        quad_f :: Quad.AbstractQuad,
                                        par :: Parameter.Parameter_MsFEM,
                                        problem_f :: Problem.AbstractBasisProblem,
                                        u_orig :: Array{Float64,1},
                                        k_time :: Int,
                                        ind_cell :: Int)
    # Abbreviations
    k = par.k
    mesh = mesh_collection.mesh
    m_f = mesh_collection.mesh_f[ind_cell]
    n_dof = m_f.n_point
    dof_f = dof_collection.dof_f[ind_cell]
    

    # ind_corner = [sort(circshift(m_f.segment[:,m_f.segment_marker.==1],1)[:])[1] ; 
    #                 sort(circshift(m_f.segment[:,m_f.segment_marker.==1],1)[:])[end] ; 
    #                 sort(circshift(m_f.segment[:,m_f.segment_marker.==2],1)[:])[end] ]
    ind_corner = [sort(m_f.segment[:,m_f.segment_marker.==1],dims=1)[1] ; 
                    sort(m_f.segment[:,m_f.segment_marker.==1],dims=1)[end] ; 
                    sort(m_f.segment[:,m_f.segment_marker.==2],dims=1)[end] ]
    # U = solution.u[mesh_collection.mesh.cell[:,ind_cell],k_time-1]
    U = u_orig[ind_corner]

    # Initialize reference basis as convected standard basis
    u0 = Problem.u_init(problem_f, m_f.point)

    # Initialize optimal basis as zeros
    uOpt = zeros(m_f.n_point,3)

    # Reconstruct edges
    reconstruct_edge_L2!(uOpt, m_f, U, u0, u_orig, 1, k[4])
    reconstruct_edge_L2!(uOpt, m_f, U, u0, u_orig, 2, k[4])
    reconstruct_edge_L2!(uOpt, m_f, U, u0, u_orig, 3, k[4])

    # Reconstruct interior
    uOpt = vec(uOpt)

    # Contraints for nodal values of basis at boundary for ALL THREE BASIS
    # FUNCTIONS!!!
    ind_con = [findall(vec(m_f.point_marker).!=0) ; 
                findall(vec(m_f.point_marker).!=0).+n_dof ; 
                findall(vec(m_f.point_marker).!=0).+2*n_dof]

    constr_val = uOpt[ind_con]
    ind_uncon = setdiff(1:(3*m_f.n_point), ind_con)

    # -------------------------------------------
    laplace_matrix = FEM.assemble_Laplace(m_f, 
                                            dof_f,
                                            ref_el_f, 
                                            quad_f)
    A = laplace_matrix' * laplace_matrix
    # A = laplace_matrix
    # -------------------------------------------

    k_reg = 0
    system_matrix = [   (U[1]^2+k_reg)*speye(n_dof)+k[1]*A      U[1]*U[2]*speye(n_dof)                      U[1]*U[3]*speye(n_dof) ; 
                        U[2]*U[1]*speye(n_dof)                  (U[2]^2+k_reg)*speye(n_dof)+k[2]*A          U[2]*U[3]*speye(n_dof) ; 
                        U[3]*U[1]*speye(n_dof)                  U[3]*U[2]*speye(n_dof)                      (U[3]^2+k_reg)*speye(n_dof)+k[3]*A    ]

    rhs = [ k_reg*u0[:,1] + U[1]*u_orig ;
            k_reg*u0[:,1] + U[2]*u_orig ;
            k_reg*u0[:,3] + U[3]*u_orig ] - system_matrix[:,ind_con]*constr_val

    uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]
    # ------------------------------------------------

    return reshape(uOpt,:,3)
end