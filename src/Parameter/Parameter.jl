module Parameter

export Parameter_FEM, Parameter_MsFEM

# -------   Type hierarchies   -------
abstract type AbstractParameter end


# --------------------------------------------------------------------------------
type Parameter_FEM <: AbstractParameter
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

    

    function Parameter_FEM(T :: Float64,
                           dt :: Float64,
                	       n_edge_per_seg :: Int64,
                           n_refinement :: Int64,
                           n_order_FEM :: Int64,
                           n_order_quad :: Int64,
                    	   time_step_method :: Int64)
        
	this = new()

	this.n_order_FEM = n_order_FEM
	this.n_edge_per_seg = n_edge_per_seg
    this.n_refinement  = n_refinement
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
type Parameter_MsFEM <: AbstractParameter
    """
    Parameter type for MsFEM.
    """
    
    n_steps :: Int64
    n_edge_per_seg :: Int64
    n_refinement :: Int64

    n_edge_per_seg_f :: Int64
    n_order_FEM_f :: Int64
    n_order_quad_f :: Int64

    T :: Float64
    dt :: Float64

    time_step_method :: Int64

    k :: Array{Float64,1}

    function Parameter_MsFEM(T :: Float64,
                             dt :: Float64,
                             n_edge_per_seg :: Int64,
                             n_refinement :: Int64,
                             n_edge_per_seg_f :: Int64,
                             n_order_FEM_f :: Int64,
                             n_order_quad_f :: Int64,
                	         time_step_method :: Int64,
                             k :: Array{Float64,1})
        
	this = new()

	this.n_edge_per_seg = n_edge_per_seg
    this.n_refinement = n_refinement
    this.n_steps = round(floor( T/dt ))
        
	this.dt = dt

	this.T = this.n_steps * dt

	this.time_step_method = time_step_method


    # Fine scale parameters
    this.n_edge_per_seg_f = n_edge_per_seg_f
    this.n_order_FEM_f = n_order_FEM_f
    this.n_order_quad_f = n_order_quad_f

       this.k = k        
	   
       return this
    end # end constructor
end # end type
# --------------------------------------------------------------------------------

end # end module
