# ---------------------------------------------------------------------------------------------
# FEM assembling
# ---------------------------------------------------------------------------------------------


function assemble_mass(solution :: Solution_FEM,
                       mesh :: Mesh.TriMesh,
                       dof :: FEM.Dof,
                       ref_el :: FEM.RefEl_Pk,
                       quad :: Quadrature.Quad_simplex,
                       par :: Parameter.Parameter_FEM,
                       problem :: Problem.Problem_top,
                       k_time :: Int64)
    
    """

    

    """

    # Assembly pattern is []
    i = d.get_dof_elem(1:d.n_elem)
    ind_test = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
    ind = vec(transpose(repmat(i, 1, size(i,2))))

    time = k_time * par.dt
    
    # Assemble element matrices in a list of size (n, n, n_elem)
    mat_local = assemble_elem_a(mesh, ref_el , quad, problem, time)

    Mat_global = sparse(ind_test, ind, vec(mat_local), dof.n_true_dof, dof.n_true_dof)
    
    
end


function assemble_elem_a(mesh :: Mesh.TriMesh,
                         ref_el :: RefEl_Pk,
                         quad :: Quadrature.Quad_simplex,
                         problem :: Problem.Problem_top,
                         time :: Float64)
    
    n = length(ref_el.n_node)
    a = zeros(n, n, mesh.n_elem)
    
    x = Mesh.map_ref_point(mesh, quad.point)
    
    velocity = problem.velocity(t, x)
    weight_elem = Grid1D.map_ref_point_der(mesh, 0.5)

    for k = 1:mesh.n_elem    
        a[:,:,k] = build_local_advection_matrix(ref_el, quad, velocity[k,:], weight_elem[k], par.conservation_form)
    end

    return a
end # end function

# ---------------------------------------------------------------------------------------------


function assemble_advection(solution :: Solution_FEM,
                            mesh :: Mesh.TriMesh,
                            dof :: FEM.Dof,
                            ref_el :: FEM.RefEl_Pk,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.Problem_top,
                            k_time :: Int64)

    """

    

    """

end


# ---------------------------------------------------------------------------------------------


function assemble_diffusion(solution :: Solution_FEM,
                            mesh :: Mesh.TriMesh,
                            dof :: FEM.Dof,
                            ref_el :: FEM.RefEl_Pk,
                            par :: Parameter.Parameter_FEM,
                            problem :: Problem.Problem_top,
                            k_time :: Int64)

    """

    

    """

end
