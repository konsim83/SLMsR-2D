module Parameter

export Parameter_FEM, Parameter_MsFEM


include("Parameter_abstract.jl")


# --------------------------------------------------------------------------------
type Parameter_FEM <: Parameter_top
    """
    Parameter type for standard FEM.
    """
    
    n_steps :: Int64
    n_edge_per_seg :: Int64
    n_order_FEM :: Int64
    n_order_quad :: Int64

    T :: Float64
    dt :: Float64

    time_step_method :: Int64

    function Parameter_FEM(T :: Float64,
                           dt :: Float64,
                	   n_edge_per_seg :: Int64,
                           n_order_FEM :: Int64,
                           n_order_quad :: Int64,
                	   time_step_method :: Int64)
        
	this = new()

	this.n_order_FEM = n_order_FEM
	this.n_edge_per_seg = n_edge_per_seg
	this.n_order_quad = n_order_quad
        this.n_steps = round(floor( T/dt ))
        
	this.dt = dt

	this.T = this.n_steps * dt

	this.time_step_method = time_step_method

	return this
    end # end constructor
end # end type
# --------------------------------------------------------------------------------



# --------------------------------------------------------------------------------
type Parameter_MsFEM <: Parameter_top
    """
    Parameter type for MsFEM.
    """
    
    n_steps :: Int64
    n_edge_per_seg :: Int64
    n_order_FEM :: Int64
    n_order_quad :: Int64

    T :: Float64
    dt :: Float64

    time_step_method :: Int64

    function Parameter_MsFEM(T :: Float64,
                             dt :: Float64,
                	     n_edge_per_seg :: Int64,
                             n_edge_per_seg_f :: Int64,
                             n_order_FEM_f :: Int64,
                             n_order_quad_f :: Int64,
                	     time_step_method :: Int64)
        
	this = new()

	this.n_edge_per_seg = n_edge_per_seg
        this.n_steps = round(floor( T/dt ))
        
	this.dt = dt

	this.T = n_steps * dt

	this.time_step_method = time_step_method


        # Fine scale parameters
        this.n_edge_per_seg_f = n_edge_per_seg_f
        this.n_order_FEM_f = n_order_FEM_f
        this.n_order_quad = n_order_quad

        
	return this
    end # end constructor
end # end type
# --------------------------------------------------------------------------------

end # end module
