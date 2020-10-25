
struct Operators{NLSIntegrator, DispFunc, FFTPlan, InvFFTPlan}
    K̂::DispFunc
    F̂::FFTPlan
    F̃̂::InvFFTPlan
    B̂::DiffEqBase.DEIntegrator
    T̂::NLSIntegrator
end # Simulation

function Operators(sim)
    # Generate FFT Plans to optimize performance
    #println("Generating FFT Plan")
    F = plan_fft!(@view sim.ψ[:, 1]) # 26 allocs
    F̃ = plan_ifft!(@view sim.ψ[:, 1]) # 34 allocs

    # Cache the kinetic factor
    @memoize function K(dx)
        println("Computing and caching K(dx = $dx)")
        if sim.α == 0
            ifftshift(cis.(dx*sim.box.ω.^2/2))
        elseif sim.α > 0
            ifftshift(cis.(dx*(sim.box.ω.^2/2 - sim.α*sim.box.ω.^3)))
        else
            throw(ArgumentError("Unknown setup α = $(sim.α) and x_order = $(sim.x_order)"))
        end
    end

    # Select algorithm
    t_algo = BS3() # Default Algorithm
    if sim.x_order == 1 && sim.α == 0
        T̂ = T₁ˢ 
    elseif sim.x_order == 2 && sim.α == 0
        T̂ = T₂ˢ
    elseif sim.x_order == 4 && sim.α == 0
        T̂ = T₄ˢ 
    elseif sim.x_order == 6 && sim.α == 0
        T̂ = T₆ˢ 
    elseif sim.x_order == 8 && sim.α == 0
        T̂ = T₈ˢ 
    elseif sim.x_order == 1 && sim.α >= 0
        T̂ = T₁ʰ 
    elseif sim.x_order == 2 && sim.α >= 0
        T̂ = T₂ʰ 
    elseif sim.x_order == 4 && sim.α >= 0
        T̂ = T₄ʰ 
        t_algo = BS3()
    else
        throw(ArgumentError("No solver available "))
    end

    # Create the integrator for the Burger term
    if sim.α > 0
        D = CenteredDifference(1, sim.t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        function MB!(du, u,p,t)
            du .= 6*sim.α*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(MB!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
    else
        function M0!(du, u,p,t)
            du .= 0*u 
        end
        prob = ODEProblem(M0!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
    end

    ops = Operators(K, F, F̃, B̂, T̂)

    return ops
end

"""
    solve!(sim::Simulation)

Solves the `Simulation` object `sim` using the techniques its attributes specify.

See also: [`init_sim`](@ref), [`NLSS.Plotter.plot_ψ`](@ref)
"""
function solve!(sim::Sim)
    #println("==========================================")
    #println("Solving cubic NLSE with the following options:")
    #print(sim)       

    # Find
    sim.ψ[:, 1] = sim.ψ₀
    # Check for pruning and calculate indices
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
    end
    
    ops = Operators(sim)

    #println("Starting evolution")
    #@showprogress 1 "Evolving in x" for i = 1:sim.box.Nₓ-1
    for i = 1:sim.box.Nₓ-1
        sim.ψ[:, i+1] .= ops.T̂(sim.ψ[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            ops.F̂*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ[j, i+1]))
            end
            F̃̂*view(sim.ψ,:,i+1)
        end # if
    end # for
    sim.solved = true

    #println("Computation Done!")
    #println("==========================================")

    return nothing
end #solve
