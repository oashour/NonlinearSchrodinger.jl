export Algorithm

"""
    solve!(sim::Simulation)

Solves the `Simulation` object `sim` using the techniques its attributes specify.

See also: [`init_sim`](@ref), [`NLSS.Plotter.plot_ψ`](@ref)
"""
function solve!(sim::Sim)
    println("==========================================")
    println("Solving cubic NLSE with the following options:")
    print(sim)       

    # Find
    sim.ψ[:, 1] = sim.ψ₀
    # Check for pruning and calculate indices
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
    end
    # Generate FFT Plans to optimize performance
    println("Generating FFT Plan")
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
    if sim.x_order == 1 && sim.α == 0
        step = T₁ˢ 
    elseif sim.x_order == 2 && sim.α == 0
        step = T₂ˢ
    elseif sim.x_order == 4 && sim.α == 0
        step = T₄ˢ 
    elseif sim.x_order == 6 && sim.α == 0
        step = T₆ˢ 
    elseif sim.x_order == 8 && sim.α == 0
        step = T₈ˢ 
    elseif sim.x_order == 1 && sim.α >= 0
        step = T₁ʰ 
        t_algo = BS3()
    elseif sim.x_order == 2 && sim.α >= 0
        step = T₂ʰ 
        algo = BS3()
    elseif sim.x_order == 4 && sim.α >= 0
        step = T₄ʰ 
        algo = BS3()
    else
        throw(ArgumentError("No solver available "))
    end

    # Create the integrator for the Burger term
    if sim.α > 0
        D = CenteredDifference(1, sim.t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        t0 = 0.0
        t1 = sim.box.dx
        function M!(du, u,p,t)
            du .= 6*sim.α*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(M!, sim.ψ[:, 1], (t0, t1))
        integrator = init(prob, algo; dt=sim.box.dx,save_everystep=false)  
    else
        integrator = 0
    end

    println("Starting evolution")
    @showprogress 1 "Evolving in x" for i = 1:sim.box.Nₓ-1
        @time sim.ψ[:, i+1] .= step(sim.ψ[:, i], K, sim.box.dx, F, F̃, integrator)
        # Pruning
        if sim.αₚ > 0 
            F*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ[j, i+1]))
            end
            F̃*view(sim.ψ,:,i+1)
        end # if
    end # for
    sim.solved = true

    println("Computation Done!")
    println("==========================================")

    return nothing
end #solve
