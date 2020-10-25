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
    if sim.α == 0
        W = ifftshift(cis.(sim.box.dx*sim.box.ω.^2/2)) # 9 allocs
    elseif sim.α > 0 && sim.x_order == 1
        W = ifftshift(cis.(sim.box.dx*(sim.box.ω.^2/2 - sim.α*sim.box.ω.^3))) # 9 allocs
    elseif sim.α > 0 && sim.x_order >= 2
        W = ifftshift(cis.(sim.box.dx/2*(sim.box.ω.^2/2 - sim.α*sim.box.ω.^3))) # 9 allocs
    else
        throw(ArgumentError("Unknown setup α = $(sim.α) and x_order = $(sim.x_order)"))
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
        integrator = init(prob, BS3(); dt=sim.box.dx,save_everystep=false)  
    else
        integrator = 0
    end

    # Select algorithm
    if sim.x_order == 1 && sim.α == 0
        step = T₁ˢ 
    elseif sim.x_order == 1 && sim.α >= 0
        step = T₁ˢʰ 
    elseif sim.x_order == 1 && sim.α == 0
        step = T₂ˢ
    elseif sim.x_order == 2 && sim.α >= 0
        step = T₁ˢʰ 
    elseif sim.x_order == 4 && sim.α == 0
        step = T₄ˢ 
    elseif sim.x_order == 6 && sim.α == 0
        step = T₆ˢ 
    elseif sim.x_order == 8 && sim.α == 0
        step = T₈ˢ 
    else
        throw(ArgumentError("No solver available "))
    end

    println("Starting evolution")
    @showprogress 1 "Evolving in x" for i = 1:sim.box.Nₓ-1
        @views sim.ψ[:, i+1] .= step(sim.ψ[:, i], W, sim.box.dx, F, F̃, integrator)
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

#######################################################################################
# Algorithms
#######################################################################################
"""
    T₁ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₁ˢʰ(ψ, W, dx, F, F̃, integrator)

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx * (-1*abs2(ψ[i]))) 
    end

    # Kinetic
    F*ψ 
    @inbounds for i in 1:length(W)
        ψ[i] *= W[i]
    end
    F̃*ψ 

    # Burger
    set_u!(integrator, ψ)
    step!(integrator)
    ψ = integrator.u

    return ψ
end #T₁ʰ
"""
    T₂ˢʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₂ˢʰ(ψ, W, dx, F, F̃, integrator)

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    # Kinetic
    F*ψ 
    @inbounds for i in 1:length(W)
        ψ[i] *= W[i]
    end
    F̃*ψ 

    # Burger
    set_u!(integrator, ψ)
    step!(integrator)
    ψ = integrator.u

    # Kinetic
    F*ψ 
    @inbounds for i in 1:length(W)
        ψ[i] *= W[i]
    end
    F̃*ψ 

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    return ψ
end #T₁ʰ
"""
    T₂ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₂ˢ(ψ, W, dx, F, F̃, integrator = 0)
    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end
    # Kinetic
    F*ψ # 0 allocs
    @inbounds for i in 1:length(W)
        ψ[i] *= W[i]
    end
    F̃*ψ # 0 allocs

    return ψ
end #T2

"""
    T₄ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T₄ˢ(ψ, W, dx, F, F̃, integrator = 0)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₂ˢ(ψ, W, ft * dx, F, F̃)
    ψ = T₂ˢ(ψ, W, bt * dx, F, F̃)
    ψ = T₂ˢ(ψ, W, ft * dx, F, F̃)

    return ψ
end # T4S

"""
    T₆ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₄ˢ`](@ref)
"""
function T₆ˢ(ψ, W, dx, F, F̃, integrator = 0)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₄ˢ(ψ, W, ft * dx, F, F̃)
    ψ = T₄ˢ(ψ, W, bt * dx, F, F̃)
    ψ = T₄ˢ(ψ, W, ft * dx, F, F̃)

    return ψ
end #T6S

"""
    T₈ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₆ˢ`](@ref)
"""
function T₈ˢ(ψ, W, dx, F, F̃, integrator = 0)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₆ˢ(ψ, W, ft * dx, F, F̃)
    ψ = T₆ˢ(ψ, W, bt * dx, F, F̃)
    ψ = T₆ˢ(ψ, W, ft * dx, F, F̃)

    return ψ
end #T8S