const N = 4000

@time begin

    A = rand(N,N)
    r = Array{Float64,1}(N)

    for j in 1:N

        s = 0.0
        for i in 1:N
            s += A[i,j]
        end
        r[j] = s
    end
end


@time begin

    A = rand(N,N)
    r = Array{Float64,1}(N)

    for i in 1:N

        s = 0.0
        for j in 1:N
            s += A[i,j]
        end
        r[i] = s
    end
end
