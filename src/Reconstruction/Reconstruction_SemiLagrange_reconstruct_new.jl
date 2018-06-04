function reconstruct_L2_conformal_new(solution :: FEM.Solution_MsFEM,
    									mesh_collection :: Mesh.TriMesh_collection,
                                        dof_collection :: FEM.AbstractDofCollection,
                                        ref_el_local :: FEM.AbstractRefEl,
                                        quad_local :: Quad.AbstractQuad,
    									par :: Parameter.Parameter_MsFEM,
    									problem :: Problem.AbstractPhysicalProblem,
    									problem_local :: Problem.AbstractBasisProblem,
    									point :: Array{Float64,2},
    									ind_cell :: Int)
    
    # Evaluate the solution at the traced back points
    uOrig = Problem.u_init(problem, point)

    uBasisOpt = reconstruct_L2_conformal_new(solution,
                                                mesh_collection,
                                                dof_collection,
                                                ref_el_local,
                                                quad_local,
                                                par,
                                                problem_local,
                                                uOrig,
                                                2,
                                                ind_cell)
    
	return uBasisOpt
 
end



# ===========================================================================================



function reconstruct_L2_conformal_new(solution :: FEM.Solution_MsFEM,
                                        mesh_collection :: Mesh.TriMesh_collection,
                                        dof_collection :: FEM.AbstractDofCollection,
                                        ref_el_local :: FEM.AbstractRefEl,
                                        quad_local :: Quad.AbstractQuad,
                                        par :: Parameter.Parameter_MsFEM,
                                        problem_local :: Problem.AbstractBasisProblem,
                                        uOrig :: Array{Float64,1},
                                        k_time :: Int,
                                        ind_cell :: Int)
    # Abbreviations
    mesh = mesh_collection.mesh
    mesh_local = mesh_collection.mesh_f[ind_cell]
    dof_local = dof_collection.dof_f[ind_cell]
    k = par.k

    # Global weights of basis
    uGlobal = solution.u[mesh.cell[:,ind_cell],k_time-1]

    # Initialize reference basis as convected standard basis
    uBasis0 = Problem.u_init(problem_local, mesh_local.point)

    # Initialize optimal basis as zeros
    uBasisOpt = zeros(mesh_local.n_point,3)

    # Reconstruct edges
    reconstruct_edge_new!(uBasisOpt, mesh_local, uGlobal, uBasis0, uOrig, 1, k[4])
    reconstruct_edge_new!(uBasisOpt, mesh_local, uGlobal, uBasis0, uOrig, 2, k[4])
    reconstruct_edge_new!(uBasisOpt, mesh_local, uGlobal, uBasis0, uOrig, 3, k[4])

    # Reconstruct interior
    uBasisOpt = vec(uBasisOpt)

    # Contraints for nodal values of basis at boundary for ALL THREE BASIS
    # FuGlobalNCTIONS!!!
    ind_con = [find(mesh_local.point_marker.!=0) ; 
                find(mesh_local.point_marker.!=0)+mesh_local.n_point ; 
                find(mesh_local.point_marker.!=0)+2*mesh_local.n_point]

    constr_val = uBasisOpt[ind_con]
    ind_uncon = setdiff(1:(3*mesh_local.n_point), ind_con)

    # System matrix and right-hand side
    n = length(uOrig)
    # system_matrix = [   (uGlobal[1]^2+k[1])*speye(n)    uGlobal[1]*uGlobal[2]*speye(n)      uGlobal[1]*uGlobal[3]*speye(n) ; 
    #                     uGlobal[2]*uGlobal[1]*speye(n)  (uGlobal[2]^2+k[2])*speye(n)        uGlobal[2]*uGlobal[3]*speye(n) ; 
    #                     uGlobal[3]*uGlobal[1]*speye(n)  uGlobal[3]*uGlobal[2]*speye(n)      (uGlobal[3]^2+k[3])*speye(n)]

    # rhs = [ k[1]*uBasis0[:,1] + uGlobal[1]*uOrig ;
    #         k[2]*uBasis0[:,2] + uGlobal[2]*uOrig ;
    #         k[3]*uBasis0[:,3] + uGlobal[3]*uOrig] - system_matrix[:,ind_con]*constr_val


    laplace_matrix = FEM.assemble_Laplace(mesh_local, 
                                            dof_local,
                                            ref_el_local, 
                                            quad_local)
    # A = laplace_matrix' * laplace_matrix
    A = laplace_matrix

    system_matrix = [   uGlobal[1]^2*speye(n)+k[1]*A    uGlobal[1]*uGlobal[2]*speye(n)      uGlobal[1]*uGlobal[3]*speye(n) ; 
                        uGlobal[2]*uGlobal[1]*speye(n)  uGlobal[2]^2*speye(n)+k[2]*A        uGlobal[2]*uGlobal[3]*speye(n) ; 
                        uGlobal[3]*uGlobal[1]*speye(n)  uGlobal[3]*uGlobal[2]*speye(n)      uGlobal[3]^2*speye(n)+k[3]*A    ]

    rhs = [ uGlobal[1]*uOrig ;
            uGlobal[2]*uOrig ;
            uGlobal[3]*uOrig] - system_matrix[:,ind_con]*constr_val

    uBasisOpt[ind_uncon] = system_matrix[ind_uncon,ind_uncon] \ rhs[ind_uncon]
    # ------------------------------------------------

    return reshape(uBasisOpt,:,3)
end