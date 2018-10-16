module Parameter

export Parameter_FEM, Parameter_MsFEM

# -------   Type hierarchies   -------
abstract type AbstractParameter end


# --------------------------------------------------------------------------------
struct Parameter_FEM <: AbstractParameter
    """
    Parameter type for standard FEM.
    """
    
    n_steps :: Int64
    n_edge_per_seg :: Int64
    n_refinement :: Int64
    n_order_FEM :: Int64
    n_order_quad :: Int64

    T :: Float64
    dt :: Float64

    time_step_method :: Int64    
end # end type


function Parameter_FEM(T :: Float64,
                       dt :: Float64,
                       n_edge_per_seg :: Int64,
                       n_refinement :: Int64,
                       n_order_FEM :: Int64,
                       n_order_quad :: Int64,
                       time_step_method :: Int64)
    
    n_steps = round(floor( T/dt ))

    Tapprox = n_steps * dt
        
    return Parameter_FEM(n_steps,
                            n_edge_per_seg,
                            n_refinement,
                            n_order_FEM,
                            n_order_quad,
                            Tapprox,
                            dt,
                            time_step_method)
end # end constructor
# --------------------------------------------------------------------------------



# --------------------------------------------------------------------------------
struct Parameter_MsFEM <: AbstractParameter
    """
    Parameter type for MsFEM.
    """
    
    n_steps :: Int64
    n_steps_f :: Int64
    n_edge_per_seg :: Int64
    n_refinement :: Int64

    n_edge_per_seg_f :: Int64
    n_order_FEM_f :: Int64
    n_order_quad_f :: Int64

    max_are_cell_f :: Float64

    T :: Float64
    dt :: Float64

    time_step_method :: Int64

    k :: Array{Float64,1}
    reconstruction_method :: Int64

    reconstruct_edge :: Bool
end # end type


function Parameter_MsFEM(T :: Float64,
                             dt :: Float64,
                             n_steps_f :: Int,
                             n_edge_per_seg :: Int64,
                             n_refinement :: Int64,
                             n_edge_per_seg_f :: Int64,
                             max_are_cell_f :: Float64,
                             n_order_FEM_f :: Int64,
                             n_order_quad_f :: Int64,
                             time_step_method :: Int64,
                             reconstruction_method :: Int64,
                             reconstruct_edge :: Bool,
                             k :: Array{Float64,1})
        
        
        n_steps = round(floor( T/dt ))
        Tapprox = n_steps * dt

        
       return Parameter_MsFEM(n_steps,
                                n_steps_f,
                                n_edge_per_seg,
                                n_refinement,
                                n_edge_per_seg_f,
                                n_order_FEM_f,
                                n_order_quad_f,
                                max_are_cell_f,
                                Tapprox,
                                dt,
                                time_step_method,
                                k,
                                reconstruction_method,
                                reconstruct_edge)
    end # end constructor
# --------------------------------------------------------------------------------

end # end module
