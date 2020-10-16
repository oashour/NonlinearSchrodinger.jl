export solve!

function solve!(sim::Simulation)
    sim.ψ[1, :] = sim.ψ₀

    println("==========================================")
    println("Solving cubic NLSE with the following options:")
    print(sim)
    @showprogress 1 "Computing..." for i = 1:sim.box.Nₓ-1
    #for i = 1:sim.box.Nₓ-1
        if sim.algorithm == "2S"
            sim.ψ[i+1, :] = T2(sim.ψ[i, :], sim.box.ω, sim.box.dx)
        elseif sim.algorithm == "4S"
            sim.ψ[i+1, :] = T4S(sim.ψ[i, :], sim.box.ω, sim.box.dx)
        elseif sim.algorithm == "6S"
            sim.ψ[i+1, :] = T6S(sim.ψ[i, :], sim.box.ω, sim.box.dx)
        elseif sim.algorithm == "8S"
            sim.ψ[i+1, :] = T8S(sim.ψ[i, :], sim.box.ω, sim.box.dx)
        else
            throw(ArgumentError("Algorithm type unknown, please check the documentation"))
        end
        # Pruning
        if sim.αₚ > 0 
            sim.ψ[i+1, :] = prune(sim.ψ[i+1, :], sim)
        end
    end #for
    sim.solved = true
    

    println("Computation Done!")
    println("==========================================")

    return nothing
end #solve

function prune(ψ, sim)
    #(naive approach)
    fft!(ψ)
    for i = 2:Int(sim.box.Nₜ/2)+1
        if (i - 1) % sim.box.n_periods != 0
           ψ[i] *= exp(-sim.αₚ*abs(ψ[i])) # + component
            ψ[sim.box.Nₜ - i + 2] *= exp(-sim.αₚ*abs(ψ[sim.box.Nₜ - i + 2])) # - component
        end
    end 
    ifft!(ψ)


    return ψ

    #   ind = []
    #   for n = 1:Int(ceil(22/2/3))
    #   for m = 1:3-1
    #   index = 3*n-m
    #   if index <  22/2
    #      append!(ind, index)
    #      append!(ind, -index)
    #   elseif index == 22/2
    #      append!(ind, -index)
    #   end
    #   end
    #   end
    #   println(sort(ind))
end


function T2(ψ, ω, dx)
    # Nonlinear
    V = -1 * abs.(ψ) .^ 2.0                      
    ψ = exp.(-im * dx / 2.0 * V) .* ψ            

    # Kinetic
    fft!(ψ)                                   
    ψ = ifftshift(exp.(-im * dx * ω .^ 2 / 2)) .* ψ  
    ifft!(ψ)                                  

    # Nonlinear
    V = -1 * abs.(ψ) .^ 2                   
    ψ = exp.(-im * dx / 2 * V) .* ψ   

    return ψ
end #T2

function T4S(ψ, ω, dx)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2(ψ, ω, ft * dx)
    ψ = T2(ψ, ω, bt * dx)
    ψ = T2(ψ, ω, ft * dx)

    return ψ
end # T4S

function T6S(ψ, ω, dx; pruning = "off", n_periods = 1)
    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T4S(ψ, ω, ft * dx)
    ψ = T4S(ψ, ω, bt * dx)
    ψ = T4S(ψ, ω, ft * dx)

    return ψ
end #T6S

function T8S(ψ, ω, dx; pruning = "off", n_periods = 1)
    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T6S(ψ, ω, ft * dx)
    ψ = T6S(ψ, ω, bt * dx)
    ψ = T6S(ψ, ω, ft * dx)

    return ψ
end #T8S