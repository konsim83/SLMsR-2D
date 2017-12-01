module Time_integrator

using FEM, Quad, Parameter, Mesh, Problem



abstract type AbstractTime_integrator end

abstract type AbstractSystem_data end
abstract type AbstractSystem_data_implEuler <: AbstractSystem_data end



# ------------------------------------------------------------------------------
# -------   ImplEuler for ADE   -------
type System_data_impEuler_ADE <: AbstractSystem_data_implEuler

    mass :: SparseMatrixCSC{Float64,Int64}
    advection :: SparseMatrixCSC{Float64,Int64}
    diffusion :: SparseMatrixCSC{Float64,Int64}
    rhs :: Array{Float64, 1}
    
    system_matrix :: SparseMatrixCSC{Float64,Int64}
    system_rhs :: Array{Float64, 1}
    
    function System_data_impEuler_ADE(dof :: FEM.AbstractDof)
        
        # Create a pattern
        i = FEM.get_dof_elem(dof, mesh, 1:dof.n_elem)
        ind = vec(i[:,[1 ; 1 ; 1 ; 2 ; 2 ; 2 ; 3 ; 3 ; 3]]')
        ind_test = vec(transpose(repmat(i, 1, size(i,2))))
        
        mass = sparse(ind_test, ind, ones(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        advection = sparse(ind_test, ind, ones(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        diffusion = sparse(ind_test, ind, ones(Float64, length(ind)), dof.n_true_dof, dof.n_true_dof)
        rhs = zeros(dof.n_true_dof)

        system_matrix = mass[dof.ind_node_non_dirichlet,dof.ind_node_non_dirichlet]
        system_rhs = FEM.map_ind_mesh2_dof(dof,rhs)
        
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
                        system_data :: System_data_impEuler_ADE,
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
