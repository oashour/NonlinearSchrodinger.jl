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
    W = ifftshift(cis.(sim.box.dx*sim.box.ω.^2/2)) # 9 allocs
    # Step through "time"
    # print("Starting evolution")
    @showprogress 1 "Evolving in x" for i = 1:sim.box.Nₓ-1
    #for i = 1:sim.box.Nₓ-1
        @views sim.ψ[:, i+1] .= sim.step(sim.ψ[:, i], W, sim.box.dx, F, F̃) # 6 allocs
        # Pruning TODO: rewrite
        if sim.αₚ > 0 
            ψ = sim.ψ[:, i+1]
            F*ψ
            @views ψ[ind_p] .*= exp.(-sim.αₚ*abs.(ψ[ind_p]))
            inv(F)*ψ
            sim.ψ[:, i+1] = ψ
        end
    end #for
    sim.solved = true

    println("Computation Done!")
    println("==========================================")

    return nothing
end #solve

#######################################################################################
# Algorithms
#######################################################################################
"""
    T₂ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₂ˢ(ψ, W, dx, F, F̃)
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

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    return ψ
end #T2

"""
    T₄ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T₄ˢ(ψ, W, dx, F, F̃)
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
function T₆ˢ(ψ, W, dx, F, F̃)

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
function T₈ˢ(ψ, W, dx, F, F̃)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₆ˢ(ψ, W, ft * dx, F, F̃)
    ψ = T₆ˢ(ψ, W, bt * dx, F, F̃)
    ψ = T₆ˢ(ψ, W, ft * dx, F, F̃)

    return ψ
end #T8S