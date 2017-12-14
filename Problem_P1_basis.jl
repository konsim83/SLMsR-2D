type BasisFun <: AbstractProblem
    
    info :: String
    type_info :: String
    
    T :: Float64

    ind_basis :: Int64
    coeff :: Array{Float64,1}

    is_transient_diffusion :: Bool
    is_transient_velocity :: Bool
    
    function BasisFun(ind_basis :: Int64, T :: Float64)

        this = new()
        
        this.info = "Evolution of basis finction with index   $ind_basis."
        this.type_info = "ADE"

        this.T = T

        this.ind_basis = ind_basis
        this.coeff = [-1.0 -1.0 1.0 ;
                      1.0 0.0 0.0 ;
                      0.0 1.0 0.0][ind_basis,:]

        is_transient_diffusion = false
        is_transient_velocity = false
        
        return this
    end # end constructor
end # end type



# --------------------------------------------------------------
# ----------------------   Functions   ----------------------
# --------------------------------------------------------------


function diffusion(problem :: BasisFun, t :: Float64, x :: Array{Float64,2})
    """

    Diffusion is represented by a positive 2-by-2 tensor.

    """
    
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2)

    out[:,1,1] = 0.01
    out[:,2,1] = 0.0
    out[:,1,2] = 0.0
    out[:,2,2] = 0.002 
    
    return out
end


function diffusion(problem :: BasisFun,  t :: Float64, x :: Array{Float64,3})
    """

    Diffusion is represented by a positive 2-by-2 tensor.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
        
    out = Array{Float64}(size(x,1), 2, 2, size(x,3))

    for i=1:size(x,3)
        out[:,:,:,i] = diffusion(problem, t, x[:,:,i])
    end
    
    return out
    
end


# --------------------------------------------------------------------


function velocity(problem :: BasisFun,  t :: Float64, x :: Array{Float64,2})
    """

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :
    
    return 0*ones(size(x))
end


function velocity(problem :: BasisFun,  t :: Float64, x :: Array{Float64,3})
    """

    Velocity is represented by a 2-vector. The solenoidal part can be
    represented by a stream function.

    """

    size(x,2)!=2 ? error(" List of vectors x must be of size nx2.") :

    #out = Array{Float64, 3}(size(x,1), 2, size(x,3))

    return 0*ones(size(x))
    
end


# --------------------------------------------------------------------

function u_init(problem :: BasisFun, x :: Array{Float64,1})
                
    length(x)!=2 ? error(" Vector x must length=2.") :
   
    out = problem.coeff' * [x ; 1] 
                
    return out
end


function u_init(problem :: BasisFun, x :: Array{Float64,2})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :
    
    out = [x ones(size(x,1))] * problem.coeff
    
    return out
end

function u_init(problem :: BasisFun, x :: Array{Float64,3})
                
    size(x,2)!=2 ? error(" List of vectors x must be of size nx2 or nx2xn_cell.") :
        
    out = Array{Float64, 3}(size(x,1), 1, size(x,3))

    for i=1:size(x,3)
        out[:,1,i] = u_init(problem, x[:,:,i])
    end
    
    return out
end
