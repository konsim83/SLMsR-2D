module Geometry

abstract type PolygonDomain end

#------------------------------------------------------------
struct Triangle <: PolygonDomain
    
    point :: Array{Float64, 2}
    segment :: Array{Int64, 2}
    segment_marker :: Array{Int64, 1}
    
    function Triangle(point :: Array{Float64,2})
        
        segment  = [1 2;
                    2 3 ;
                    3 1]
        
        segment_marker = [1 ; 2 ; 3]

        return new(point, segment, segment_marker)
    end # end constructor    
end
#------------------------------------------------------------

#------------------------------------------------------------
struct UnitSquare <: PolygonDomain
    
    point :: Array{Float64, 2}
    segment :: Array{Int64, 2}
    segment_marker :: Array{Int64, 1}
    
    function UnitSquare(segment_marker = [1 ; 2 ; 3 ; 4] :: Array{Int64,1})

        size(segment_marker)!=(4,1) ? error("Segment marker list must be a 4-Vector.") :
        
        point = [0.0 0.0 ;
                 1.0 0.0 ;
                 1.0 1.0 ;
                 0.0 1.0]
        
        segment  = [1 2;
                    2 3 ;
                    3 4 ;
                    4 1]
        
        return new(point, segment, segment_marker)        
    end # end constructor    
end # end type
#------------------------------------------------------------


#------------------------------------------------------------
type UnitSimplex <: PolygonDomain

    point :: Array{Float64, 2}
    segment :: Array{Int64, 2}
    segment_marker :: Array{Int64, 1}
    
    function UnitSimplex()

        point = [0.0 0.0 ;
                 1.0 0.0 ;
                 0.0 1.0]
        
        segment  = [1 2;
                    2 3 ;
                    3 1]
        
        segment_marker = [1 ; 2 ; 3]

        return new(point, segment, segment_marker)
    end # end constructor
end # end type
#------------------------------------------------------------


# -----------------------------------------------------------
# -------   Functions   ----------------------------------
# -----------------------------------------------------------

function compute_P1_basis_coeff(tri :: Triangle)
    
    # columns are coefficients
    coeff = [tri.point ones(3)] \ eye(3)
    
    return coeff
end

# -----------------------------------------------------------
# -----------------------------------------------------------
# -----------------------------------------------------------

end # end module
