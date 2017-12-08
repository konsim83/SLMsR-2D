module Time_integrator

using FEM, Quad, Parameter, Mesh, Problem



abstract type AbstractTime_integrator end

abstract type AbstractSystem_data end
abstract type AbstractSystem_data_implEuler <: AbstractSystem_data end



# ------------------------------------------------------------------------------
# -------   ImplEuler for ADE   -------
type System_data_implEuler_ADE <: AbstractSystem_data_implEuler

    mass :: SparseMatrixCSC{Float64,Int64}
    advection :: SparseMatrixCSC{Float64,Int64}
    diffusion :: SparseMatrixCSC{Float64,Int64}
    rhs :: Array{Float64, 1}
    
    system_matrix :: SparseMatrixCSC{Float64,Int64}
    system_rhs :: Array{Float64, 1}

    # holds solution of system in dof-indices. Note that Dirichlet
    # nodes can occur with periodic boundaries.
    u_temp :: Array{Float64, 1}

    #weight_elem :: Array{Float64,2}
    #DPhi :: Array{Float64,4}
    
    function System_data_implEuler_ADE(dof :: FEM.AbstractDof, mesh :: Mesh.TriMesh, problem :: Problem.AbstractProblem)
        
        # Create a pattern
        i = FEM.get_dof_elem(dof, mesh, 1:dof.n_elem)
        ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
        ind_test = vec(transpose(repmat(i, 1, size(i,2))))

        # Allocate memory for sparse matrix
        mass = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        advection = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        diffusion = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        rhs = zeros(dof.n_true_dof)

        system_matrix = mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
        system_rhs = rhs[dof.ind_node_non_dirichlet]

        #weight_elem = weight_elem = Mesh.map_ref_point_grad_det(m, q.point, 1:d.n_elem)
        
        u_temp = FEM.map_ind_mesh2dof(dof, Problem.u_init(problem, mesh.point))
        
        return new(mass, advection, diffusion, rhs, system_matrix, system_rhs, u_temp)
    end
    
end
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# -------   Implicit Euler   -------
type ImplEuler{T <: AbstractSystem_data_implEuler} <: AbstractTime_integrator

    system_data :: T
    
    function ImplEuler{T}(dof :: FEM.AbstractDof, mesh :: Mesh.TriMesh, problem :: Problem.AbstractProblem) where {T <: AbstractSystem_data_implEuler}
        
	# Reserve Memory for System data only
        if problem.type_info=="ADE"
	    system_data = System_data_implEuler_ADE(dof, mesh, problem)
        else
            error("Problem type not implemented yet.")
        end
            
        return new(system_data)
    end
end
# ------------------------------------------------------------------------------







# ------------------------------------------------------------------------------
function make_step!(solution :: FEM.AbstractSolution,
                    time_int :: ImplEuler,
                    mesh :: Mesh.TriMesh,
                    dof :: FEM.AbstractDof,
                    ref_el :: FEM.AbstractRefEl,
                    quad :: Quad.AbstractQuad,
                    par :: Parameter.AbstractParameter,
                    problem :: Problem.AbstractProblem,
                    k_time :: Int64)
    

    """

    Accept as input an Implicit Euler system in form of a type
    and do time step from k_time -> k_time+1.

    """


    # ---------------------------------------
    # Manipulate the system data locally
    update_system!(solution,
                   time_int.system_data,
                   mesh,
                   dof,
                   ref_el,
                   quad,
                   par,
                   problem,
                   k_time+1)
    # ---------------------------------------

    
    # Solve the system and map to mesh variables
    time_int.system_data.u_temp[dof.ind_node_non_dirichlet] = time_int.system_data.system_matrix \ time_int.system_data.system_rhs
    solution.u[:,k_time+1] = FEM.map_ind_dof2mesh(dof, time_int.system_data.u_temp)
    
    #solution.u[dof.ind_node_non_dirichlet,k_time+1] = FEM.map_ind_dof2mesh( dof,  )
    #solution.u[dof.ind_node_dirichlet,k_time+1] = solution.u[dof.ind_node_dirichlet, k_time]    
    return nothing
end
# ----------------------------------


function update_system!(solution :: FEM.AbstractSolution,
                        system_data :: System_data_implEuler_ADE,
                        mesh :: Mesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.AbstractRefEl,
                        quad :: Quad.AbstractQuad,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractProblem,
                        k_time :: Int64)

    #=
    # ----   This is the fast version   ---
    FEM.assemble_mass!(system_data.mass,
                       solution,
                       mesh,
                       dof,
                       ref_el,
                       quad,
                       par,
                       problem,
                       k_time)
    
    FEM.assemble_advection!(system_data.advection,
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            quad,
                            par,
                            problem,
                            k_time)
    
    FEM.assemble_diffusion!(system_data.diffusion,
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            quad,
                            par,
                            problem,
                            k_time)
    # ----   This is the fast version   ---
    =#


    
    # ----   This is the slow version   ---
    system_data.mass = FEM.assemble_mass(solution,
                                         mesh,
                                         dof,
                                         ref_el,
                                         quad,
                                         par,
                                         problem,
                                         k_time)
    
    system_data.advection = FEM.assemble_advection(solution,
                                                   mesh,
                                                   dof,
                                                   ref_el,
                                                   quad,
                                                   par,
                                                   problem,
                                                   k_time)
    
    system_data.diffusion = FEM.assemble_diffusion(
                                                   solution,
                                                   mesh,
                                                   dof,
                                                   ref_el,
                                                   quad,
                                                   par,
                                                   problem,
                                                   k_time)
    # ----   This is the slow version   ---


    
    if dof.n_node_neumann > 0
        error("Neumann boundary integral not implemented yet.")
    end

    system_data.system_matrix = ( (system_data.mass - par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] )
    
    if dof.is_periodic
        system_data.system_rhs = system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] * (FEM.map_ind_mesh2dof(dof, solution.u[:,k_time-1])[dof.ind_node_non_dirichlet])
    else
        
        system_data.system_rhs = (   (system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
                                      * solution.u[dof.ind_node_non_dirichlet,k_time-1] -
                                      (system_data.mass -
                                       par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_dirichlet]
                                      * solution.u[dof.ind_node_dirichlet,k_time-1]) )
    end
    
end


end # end module




#=


module Time_integrator

using FEM, Quad, Parameter, Mesh, Problem



abstract type AbstractTime_integrator end

abstract type AbstractSystem_data end
abstract type AbstractSystem_data_implEuler <: AbstractSystem_data end



# ------------------------------------------------------------------------------
# -------   ImplEuler for ADE   -------
type System_data_implEuler_ADE <: AbstractSystem_data_implEuler

    mass :: SparseMatrixCSC{Float64,Int64}
    advection :: SparseMatrixCSC{Float64,Int64}
    diffusion :: SparseMatrixCSC{Float64,Int64}
    rhs :: Array{Float64, 1}
    
    system_matrix :: SparseMatrixCSC{Float64,Int64}
    system_rhs :: Array{Float64, 1}
    
    function System_data_implEuler_ADE(dof :: FEM.AbstractDof)
        
        # Create a pattern
        i = FEM.get_dof_elem(dof, mesh, 1:dof.n_elem)
        ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
        ind_test = vec(transpose(repmat(i, 1, size(i,2))))
        
        mass = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        advection = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        diffusion = sparse(ind_test, ind, zeros(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        rhs = zeros(dof.n_true_dof)

        system_matrix = mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
        system_rhs = rhs[dof.ind_node_non_dirichlet]

        return new(mass, advection, diffusion, rhs, system_matrix, system_rhs)
    end
    
end
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# -------   Implicit Euler   -------
type ImplEuler{T <: AbstractSystem_data_implEuler} <: AbstractTime_integrator

    system_data :: T
    
    function ImplEuler{T}(dof :: FEM.AbstractDof, type_info :: String) where {T <: AbstractSystem_data_implEuler}
        
	# Reserve Memory for System data only
        if type_info=="ADE"
	    system_data = System_data_implEuler_ADE(dof)
        else
            error("Problem type not implemented yet.")
        end
            
        return new(system_data)
    end
end
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
function make_step!(solution :: FEM.AbstractSolution,
                    time_int :: ImplEuler,
                    mesh :: Mesh.TriMesh,
                    dof :: FEM.AbstractDof,
                    ref_el :: FEM.AbstractRefEl,
                    quad :: Quad.AbstractQuad,
                    par :: Parameter.AbstractParameter,
                    problem :: Problem.AbstractProblem,
                    k_time :: Int64)
    

    """

    Accept as input an Implicit Euler system in form of a type
    and do time step from k_time -> k_time+1.

    """


    # ---------------------------------------
    # Manipulate the system data locally
    update_system!(time_int.system_data,
                   solution,
                   mesh,
                   dof,
                   ref_el,
                   quad,
                   par,
                   problem,
                   k_time+1)    
    # ---------------------------------------
    
    
    # Solve the system and map to mesh variables
    solution.u[dof.ind_node_non_dirichlet,k_time+1] = dof.map_ind_dof2mesh( time_int.system_data.system_matrix \ time_int.system_data.system_rhs )
    solution.u[dof.ind_node_dirichlet,k_time+1] = solution.u[dof.ind_node_non_dirichlet, k_time]
    
    return nothing
end
# ----------------------------------


function update_system!(solution :: FEM.AbstractSolution,
                        system_data :: System_data_implEuler_ADE,
                        mesh :: Mesh.TriMesh,
                        dof :: FEM.AbstractDof,
                        ref_el :: FEM.AbstractRefEl,
                        quad :: Quad.AbstractQuad,
                        par :: Parameter.AbstractParameter,
                        problem :: Problem.AbstractProblem,
                        k_time :: Int64)

    FEM.assemble_mass!(system_data.mass,
                       solution,
                       mesh,
                       dof,
                       ref_el,
                       quad,
                       par,
                       problem,
                       k_time)
    
    FEM.assemble_advection!(system_data.advection,
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            quad,
                            par,
                            problem,
                            k_time)
    
    FEM.assemble_diffusion!(system_data.diffusion,
                            solution,
                            mesh,
                            dof,
                            ref_el,
                            quad,
                            par,
                            problem,
                            k_time)
    
    if dof.n_node_neumann > 0
        error("Neumann boundary integral not implemented yet.")
    end
    
    system_data.system_matrix = ( (system_data.mass
                                            -
                                          par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
                                          )
    
    if dof.is_periodic
            system_data.system_rhs = system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet] * dof.map_ind_mesh2_dof(solution.u[:,k_time-1])
    else
        
        system_data.system_rhs = (   (system_data.mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
                                      * solution.u[dof.ind_node_non_dirichlet,k_time-1] -
                                      (system_data.mass -
                                       par.dt*(system_data.diffusion-system_data.advection))[dof.ind_node_non_dirichlet,dof.ind_node_dirichlet]
                                      * solution.u[dof.ind_node_dirichlet,k_time-1]) )
    end
    
end


end # end module



=#
