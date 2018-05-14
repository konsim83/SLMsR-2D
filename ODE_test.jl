using DifferentialEquations, PyPlot

N = 10000
x = rand(-20.0:0.001:2.0,2,N)
t = (0.0, 0.01)

# Version for n systems

v = function(x, p, t)

	# du = -1.0*ones(size(x))

	dx = [0.0 1.0 ; -1.0 0.0]*x

	# du = [cos(u[2]-u[1]) ; 
	#      	sin(u[1]-u[2])]

end

function velocity(x, p, t)

    size(x,1)!=2 ? error("List of vectors x must be of size 2-by-n.") :

    out = [[sin(x[2,i]) ; sin(x[1,i])] for i=1:size(x,2)]
    
    return hcat(out...)
end


u0 = [1.0 ; 2.0]
prob = ODEProblem(velocity, x, t)
# sol = solve(prob,alg_hints=[:stiff])
sol = solve(prob,ode45())



# U = vcat(sol.u...)
# for i=1:maximum([N, 400])
# 	scatter(U[1,i], U[2,i], s=50); plot(U[1:2:end,i], U[2:2:end,i])
# end