module FiniteDiff

function central!(phi_t :: Array{Float64,2},
                  phi :: Array{Float64,2},
                  dt :: Float64)
    
    phi_t[:,2:end-1] = (phi[:,3:end] - phi[:,1:end-2]) / (2*dt)
    phi_t[:,1] = (phi[:,2] - phi[:,1]) / dt
    phi_t[:,end] = (phi[:,end] - phi[:,end-1]) / dt
    
end

end
