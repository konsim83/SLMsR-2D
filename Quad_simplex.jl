type Quad_simplex

    order :: Int64

    weight :: Array{Float64, 1}
    point:: Array{Float64, 2}

    function Quad_simplex(number :: Int64)
        this = new()

        if number==1

            this.weight = [0.5]
	    this.point = [1/3 ; 1/3]
	    this.order = 1
            
        elseif number==2

            this.weight = [1/6 ; 1/6 ; 1/6]
	    this.point = [0.5 0 ; 0.5 0.5 ; 0 0.5]
	    this.order = 2
            
        elseif number==3

            this.weight = [1/6 ; 1/6 ; 1/6]
	    this.point = [1/6 1/6 ; 1/6 4/6 ; 4/6 1/6]
	    this.order = 2
            
        elseif number==4

	    this.weight = [-9/32 ; 25/96 ; 25/96 ; 25/96]
	    this.point = [1/3 1/3 ; 0.2 0.2 ; 0.6 0.2 ; 0.2 0.6]
	    this.order = 3
            
        elseif number==5

            A = 0.109951743655322/2
	    B = 0.223381589678011/2
	    r = 0.091576213509771
	    s = 0.445948490915965
	    t = 0.816847572980459
	    u = 0.108103018168070
            
	    this.weight = [A ; A ; A ; B ; B ; B]
	    this.point = [r r ; t r ; r t ; s s ; u s ; s u]
	    this.order = 4
            
        elseif number==6
            
            A = 9/80
	    B = (155-sqrt(15))/2400
	    C = (155+sqrt(15))/2400
	    t = 1/3
	    r = (6-sqrt(15))/21
	    s = (9+2*sqrt(15))/21
	    u = (6+sqrt(15))/21
	    v = (9-2*sqrt(15))/21

	    this.weight = [A ; B ; B ; B ; C ; C ; C]
	    this.point = [t t ; r r ; r s ; s r ; u u ; u v ; v u]
	    this.order = 5

        else
            error("Quadrature rule not implemented")
        end # end if

        return this
    end # end constructor
end # end immutable type
