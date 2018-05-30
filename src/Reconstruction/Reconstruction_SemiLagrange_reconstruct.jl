function reconstruct_L2_initial(solution :: FEM.Solution_MsFEM,
									mesh_collection :: Mesh.TriMesh_collection,
									par :: Parameter.Parameter_MsFEM,
									problem :: Problem.AbstractPhysicalProblem,
									problem_f :: Problem.AbstractBasisProblem,
									point :: Array{Float64,2},
									ind_cell :: Int)
    
    # Evaluate the solution at the traced back points
    u_orig = Problem.u_init(problem, point)

    uOpt = reconstruct_L2(solution,
                            mesh_collection,
                            par,
                            problem_f,
                            u_orig,
                            2,
                            ind_cell)
    
	return uOpt
 
end



# ===========================================================================================



function reconstruct_L2(solution :: FEM.Solution_MsFEM,
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
    U = solution.u[mesh_collection.mesh.cell[:,ind_cell],k_time-1]
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
    reconstruct_edge!(uOpt, m_f, U, u0, u, 1, par.k[4])
    reconstruct_edge!(uOpt, m_f, U, u0, u, 2, par.k[4])
    reconstruct_edge!(uOpt, m_f, U, u0, u, 3, par.k[4])
    # ------------------------------------------------


    # ------------------------------------------------
    # Reconstruct interior
    uOpt = vec(uOpt)

    # Contraints for nodal values of basis at boundary for ALL THREE BASIS
    # FUNCTIONS!!!
    ind_con = [find(m_f.point_marker.!=0) ; find(m_f.point_marker.!=0)+n_dof ; find(m_f.point_marker.!=0)+2*n_dof]

    constr_val = uOpt[ind_con]
    ind_uncon = setdiff(1:(3*n_dof), ind_con)

    n = length(u)
    system_matrix = [(u1^2+k[1])*speye(n) u1*u2*speye(n) u1*u3*speye(n) ; 
                        u2*u1*speye(n) (u2^2+k[2])*speye(n) u2*u3*speye(n) ; 
                        u3*u1*speye(n) u3*u2*speye(n) (u3^2+k[3])*speye(n)]
    rhs = [k[1]*u0[:,1] + u1*u; 
            k[2]*u0[:,2] + u2*u ; 
            k[3]*u0[:,3] + u3*u]  - system_matrix[:,ind_con]*constr_val

    uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]
    # ------------------------------------------------

    return reshape(uOpt,:,3)
end




# -------------------------------------------------

# Old state before reimplementation

# function reconstruct_L2_initial(solution :: FEM.Solution_MsFEM,
#                                     mesh_collection :: Mesh.TriMesh_collection,
#                                     par :: Parameter.Parameter_MsFEM,
#                                     problem :: Problem.AbstractPhysicalProblem,
#                                     problem_f :: Problem.AbstractBasisProblem,
#                                     point :: Array{Float64,2},
#                                     ind_cell :: Int)
    
#     # ------------------------------------------------
#     # Evaluate the solution at the traced back points
#     u = Problem.u_init(problem, point)

#     U = solution.u[mesh_collection.mesh.cell[:,ind_cell],1]
#     u1 = U[1]
#     u2 = U[2]
#     u3 = U[3]
#     # ------------------------------------------------

#     u0 = Problem.u_init(problem_f, mesh_collection.mesh_f[ind_cell].point)

#     k = par.k
#     m_f = mesh_collection.mesh_f[ind_cell]
#     n_dof = m_f.n_point
#     uOpt = zeros(n_dof,3)

#     # ------------------------------------------------
#     # Reconstruct edges
    

#     # Indices of points on certain boundary edge (indices in terms
#     # of original local mesh)
#     ind_edge = sort(unique(m_f.segment[:,m_f.segment_marker.==1])) # 2-by-n matrix
#     n = length(ind_edge)

#     basis_lr = zeros(2*n)

#     u_before = u[ind_edge]
#     basis_left = u0[ind_edge,1]
#     basis_right = u0[ind_edge,2]

#     ind_con = [1 ; 2 ; 1+n ; 2+n]
#     val_con = [1. ; 0. ; 0. ; 1.]
#     ind_uncon = setdiff(1:(2*n), ind_con)

#     a_1 = u1
#     a_2 = u2

#     system_matrix = [ (a_1^2 + k[1])*speye(n)   (a_1*a_2)*speye(n) ; 
#                         (a_1*a_2)*speye(n)   (a_2^2 + k[1])*speye(n) ]

#     rhs = ( [k[1]*basis_left + a_1*u_before ;
#             k[1]*basis_right + a_2*u_before]
#             - system_matrix[:,ind_con]*val_con )
    
#     basis_lr[ind_con] = val_con
#     basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]

#     uOpt[ind_edge,1] = basis_lr[1:n]
#     uOpt[ind_edge,2] = basis_lr[(n+1):end]

#     ########################################

#     ind_edge = sort(unique(m_f.segment[:,m_f.segment_marker.==2])) # 2-by-n matrix
#     n = length(ind_edge)

#     basis_lr = zeros(2*n)

#     u_before = u[ind_edge]
#     basis_left = u0[ind_edge,2]
#     basis_right = u0[ind_edge,3]

#     ind_con = [1 ; 2 ; 1+n ; 2+n]
#     val_con = [1. ; 0. ; 0. ; 1.]
#     ind_uncon = setdiff(1:(2*n), ind_con)

#     a_1 = u2
#     a_2 = u3

#     system_matrix = [ (a_1^2 + k[1])*speye(n)   (a_1*a_2)*speye(n) ; 
#                         (a_1*a_2)*speye(n)   (a_2^2 + k[1])*speye(n) ]

#     rhs = ( [k[1]*basis_left + a_1*u_before ;
#             k[1]*basis_right + a_2*u_before]
#             - system_matrix[:,ind_con]*val_con )
    
#     basis_lr[ind_con] = val_con
#     basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]

#     uOpt[ind_edge,2] = basis_lr[1:n]
#     uOpt[ind_edge,3] = basis_lr[(n+1):end]

#     ########################################

#     ind_edge = sort(unique(m_f.segment[:,m_f.segment_marker.==3])) # 2-by-n matrix
#     n = length(ind_edge)

#     basis_lr = zeros(2*n)

#     u_before = u[ind_edge]
#     basis_left = u0[ind_edge,3]
#     basis_right = u0[ind_edge,1]

#     ind_con = [1 ; 2 ; 1+n ; 2+n]
#     val_con = [1. ; 0. ; 0. ; 1.]
#     ind_uncon = setdiff(1:(2*n), ind_con)

#     a_1 = u3
#     a_2 = u1

#     system_matrix = [ (a_1^2 + k[1])*speye(n)   (a_1*a_2)*speye(n) ; 
#                         (a_1*a_2)*speye(n)   (a_2^2 + k[1])*speye(n) ]

#     rhs = ( [k[1]*basis_left + a_1*u_before ;
#             k[1]*basis_right + a_2*u_before]
#             - system_matrix[:,ind_con]*val_con )
    
#     basis_lr[ind_con] = val_con
#     basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]

#     uOpt[ind_edge,3] = basis_lr[1:n]
#     uOpt[ind_edge,1] = basis_lr[(n+1):end]

#     # ------------------------------------------------


#     # ------------------------------------------------
#     # Reconstruct interior

#     uOpt = vec(uOpt)

#     # Contraints for nodal values of basis
#     ind_con = sort(unique(m_f.segment[:,[find(m_f.segment_marker.==1) find(m_f.segment_marker.==2) find(m_f.segment_marker.==3)]])) # 2-by-n matrix
#     constr_val = uOpt[ind_con]
#     ind_uncon = setdiff(1:(3*n_dof), ind_con)

#     n = length(u)
#     system_matrix = [(u1^2+k[1])*speye(n) u1*u2*speye(n) u1*u3*speye(n) ; 
#                         u2*u1*speye(n) (u2^2+k[2])*speye(n) u2*u3*speye(n) ; 
#                         u3*u1*speye(n) u3*u2*speye(n) (u3^2+k[3])*speye(n)]
#     rhs = [k[1]*u0[:,1] + u1*u; 
#             k[2]*u0[:,2] + u2*u ; 
#             k[3]*u0[:,3] + u3*u]  - system_matrix[:,ind_con]*constr_val

#     uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]
#     # ------------------------------------------------

#     return reshape(uOpt,:,3)
 

#     return reshape(uOpt,:,3)
# end



# # ===========================================================================================



# function reconstruct_L2(solution :: FEM.Solution_MsFEM,
#                                     mesh_collection :: Mesh.TriMesh_collection,
#                                     par :: Parameter.Parameter_MsFEM,
#                                     problem_f :: Problem.AbstractBasisProblem,
#                                     point_orig :: Array{Float64,2},
#                                     u_orig :: Array{Float64,1},
#                                     k_time :: Int,
#                                     ind_cell :: Int)

#     # ------------------------------------------------
#     # Evaluate the solution at the traced back points
#     # u = PostProcess.evaluate(solution, mesh_collection, point_orig, k_time-1)
#     u = u_orig
#     U = solution.u[mesh_collection.mesh.cell[:,ind_cell],k_time-1]
#     u1 = U[1]
#     u2 = U[2]
#     u3 = U[3]
#     # ------------------------------------------------

#     u0 = Problem.u_init(problem_f, mesh_collection.mesh_f[ind_cell].point)

#     k = par.k
#     m_f = mesh_collection.mesh_f[ind_cell]
#     n_dof = m_f.n_point
#     uOpt = zeros(n_dof,3)

#     # ------------------------------------------------
#     # Reconstruct edges
    

#     # Indices of points on certain boundary edge (indices in terms
#     # of original local mesh)
#     ind_edge = sort(unique(m_f.segment[:,m_f.segment_marker.==1])) # 2-by-n matrix
#     n = length(ind_edge)

#     basis_lr = zeros(2*n)

#     u_before = u[ind_edge]
#     basis_left = u0[ind_edge,1]
#     basis_right = u0[ind_edge,2]

#     ind_con = [1 ; 2 ; 1+n ; 2+n]
#     val_con = [1. ; 0. ; 0. ; 1.]
#     ind_uncon = setdiff(1:(2*n), ind_con)

#     a_1 = u1
#     a_2 = u2

#     system_matrix = [ (a_1^2 + k[1])*speye(n)   (a_1*a_2)*speye(n) ; 
#                         (a_1*a_2)*speye(n)   (a_2^2 + k[1])*speye(n) ]

#     rhs = ( [k[1]*basis_left + a_1*u_before ;
#             k[1]*basis_right + a_2*u_before]
#             - system_matrix[:,ind_con]*val_con )
    
#     basis_lr[ind_con] = val_con
#     basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]

#     uOpt[ind_edge,1] = basis_lr[1:n]
#     uOpt[ind_edge,2] = basis_lr[(n+1):end]

#     ########################################

#     ind_edge = sort(unique(m_f.segment[:,m_f.segment_marker.==2])) # 2-by-n matrix
#     n = length(ind_edge)

#     basis_lr = zeros(2*n)

#     u_before = u[ind_edge]
#     basis_left = u0[ind_edge,2]
#     basis_right = u0[ind_edge,3]

#     ind_con = [1 ; 2 ; 1+n ; 2+n]
#     val_con = [1. ; 0. ; 0. ; 1.]
#     ind_uncon = setdiff(1:(2*n), ind_con)

#     a_1 = u2
#     a_2 = u3

#     system_matrix = [ (a_1^2 + k[1])*speye(n)   (a_1*a_2)*speye(n) ; 
#                         (a_1*a_2)*speye(n)   (a_2^2 + k[1])*speye(n) ]

#     rhs = ( [k[1]*basis_left + a_1*u_before ;
#             k[1]*basis_right + a_2*u_before]
#             - system_matrix[:,ind_con]*val_con )
    
#     basis_lr[ind_con] = val_con
#     basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]

#     uOpt[ind_edge,2] = basis_lr[1:n]
#     uOpt[ind_edge,3] = basis_lr[(n+1):end]

#     ########################################

#     ind_edge = sort(unique(m_f.segment[:,m_f.segment_marker.==3])) # 2-by-n matrix
#     n = length(ind_edge)

#     basis_lr = zeros(2*n)

#     u_before = u[ind_edge]
#     basis_left = u0[ind_edge,3]
#     basis_right = u0[ind_edge,1]

#     ind_con = [1 ; 2 ; 1+n ; 2+n]
#     val_con = [1. ; 0. ; 0. ; 1.]
#     ind_uncon = setdiff(1:(2*n), ind_con)

#     a_1 = u3
#     a_2 = u1

#     system_matrix = [ (a_1^2 + k[1])*speye(n)   (a_1*a_2)*speye(n) ; 
#                         (a_1*a_2)*speye(n)   (a_2^2 + k[1])*speye(n) ]

#     rhs = ( [k[1]*basis_left + a_1*u_before ;
#             k[1]*basis_right + a_2*u_before]
#             - system_matrix[:,ind_con]*val_con )
    
#     basis_lr[ind_con] = val_con
#     basis_lr[ind_uncon] = system_matrix[ind_uncon, ind_uncon] \ rhs[ind_uncon]

#     uOpt[ind_edge,3] = basis_lr[1:n]
#     uOpt[ind_edge,1] = basis_lr[(n+1):end]

#     # ------------------------------------------------


#     # ------------------------------------------------
#     # Reconstruct interior

#     uOpt = vec(uOpt)

#     # Contraints for nodal values of basis
#     ind_con = sort(unique(m_f.segment[:,[find(m_f.segment_marker.==1) find(m_f.segment_marker.==2) find(m_f.segment_marker.==3)]])) # 2-by-n matrix
#     constr_val = uOpt[ind_con]
#     ind_uncon = setdiff(1:(3*n_dof), ind_con)

#     n = length(u)
#     system_matrix = [(u1^2+k[1])*speye(n) u1*u2*speye(n) u1*u3*speye(n) ; 
#                         u2*u1*speye(n) (u2^2+k[2])*speye(n) u2*u3*speye(n) ; 
#                         u3*u1*speye(n) u3*u2*speye(n) (u3^2+k[3])*speye(n)]
#     rhs = [k[1]*u0[:,1] + u1*u; 
#             k[2]*u0[:,2] + u2*u ; 
#             k[3]*u0[:,3] + u3*u]  - system_matrix[:,ind_con]*constr_val

#     uOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]
#     # ------------------------------------------------

#     return reshape(uOpt,:,3)
# end