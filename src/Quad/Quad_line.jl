struct Quad_line <: AbstractQuad
    # Gauss-Jacobi quadrature. For standard Gauss quadrature set
    # alpha=beta=0
    
    order :: Int64

    weight :: Array{Float64, 1}
    point :: Array{Float64, 1}

    n_point :: Int

    function Quad_line(alpha :: Float64,
                        beta :: Float64,
                        N :: Int)
            
        if N==0
            x = (alpha-beta) / (alpha+beta+2.0)
            w = 2
            return x, w
        end

        J = zeros(N+1)

        h1 = 2*(0:N) .+ alpha .+ beta

        var1 = vec(2. / (h1[1:N].+2))
        var2 = sqrt.( (1:N) .* ((1:N).+alpha.+beta) .* ((1:N).+alpha) .* ((1:N).+beta) ./ (h1[1:N].+1) ./ (h1[1:N].+3) )
        
        J1 = diagm(0=>-0.5*(alpha^2-beta^2) ./ (h1.+2) ./ h1)
        J2 = diagm(1=>var1.*var2)
        J = J1 + J2

        if (alpha + beta) < 10*eps()
            J[1,1] = 0
        end

        J = J + J'

        x, V = eigen(J)
        w = (V[1,:]').^2 * 2^(alpha + beta + 1)/(alpha + beta + 1) * gamma(alpha + 1)*gamma(beta + 1)/gamma(alpha + beta  + 1)
        order = 2*length(w) - 1

        return new(order, 0.5*collect(vec(w)), 0.5*(collect(vec(x)).+1), length(w))
    end # end function    
end # end immutable type
