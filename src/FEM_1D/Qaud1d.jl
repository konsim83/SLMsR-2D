struct Quad1D
    point
    weight
    order

    function Quad_jacobi_gauss(alpha, beta, N)
        
        if N==0
            x = (alpha-beta) / (alpha+beta+2.0)
            w = 2
            return x, w
        end

        J = zeros(N+1)

        h1 = 2*(0:N) + alpha + beta

        var1 = 2./(h1[1:N]+2) 
        var2 = sqrt.( (1:N).*   ((1:N)+alpha+beta).*   ((1:N)+alpha).*   ((1:N)+beta)./   (h1[1:N]+1)./   (h1[1:N]+3) )
        J = diagm(-0.5*(alpha^2-beta^2)./(h1+2)./h1) +
        diagm(   var1.*var2, 1)

        if (alpha + beta) < 10*eps()
            J[1,1] = 0
        end

        J = J + J'

        x, V = eig(J)
        w = (V[1,:]').^2 * 2^(alpha + beta + 1)/(alpha + beta + 1) * gamma(alpha + 1)*gamma(beta + 1)/gamma(alpha + beta  + 1)
        order = 2*length(w) - 1

        return new(0.5*(collect(vec(x))+1), 0.5*collect(vec(w)), order)
    end # end function
end # end type