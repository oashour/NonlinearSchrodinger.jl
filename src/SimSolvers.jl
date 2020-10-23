export Algorithm

"""
    solve!(sim::Simulation)

Solves the `Simulation` object `sim` using the techniques its attributes specify.

See also: [`init_sim`](@ref), [`NLSS.Plotter.plot_ψ`](@ref)
"""
function solve!(sim::Sim)
    #println("==========================================")
    #println("Solving cubic NLSE with the following options:")
    # Copy in x = 0 array
    sim.ψ[:, 1] = sim.ψ₀
    # Check for pruning and calculate indices
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
    end
    # Generate FFT Plans to optimize performance
    #println("Generating FFT Plan")
    F = plan_fft!(sim.ψ[:, 1]) # Plan
    F̃ = plan_ifft!(sim.ψ[:, 1]) # Plan

    # Print info about simulation
    #print(sim)

    # Cache the kinetic factor
    W = ifftshift(cis.(sim.box.dx * sim.box.ω .^ 2 / 2))
    # Step through time
    #@showprogress 1 "Computing..." for i = 1:sim.box.Nₓ-1
    for i = 1:sim.box.Nₓ-1
        sim.ψ[:, i+1] .= sim.step(sim.ψ[:, i], W, sim.box.dx, F, F̃)
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

    #println("Computation Done!")
    #println("==========================================")

    return nothing
end #solve

#######################################################################################
# Algorithms
#######################################################################################
"""
    T2(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T2(ψ, W, dx, F, F̃)
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

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    return ψ
end #T2

"""
    T4S(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T4S(ψ, W, dx, F, F̃)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2(ψ, W, ft * dx, F, F̃)
    ψ = T2(ψ, W, bt * dx, F, F̃)
    ψ = T2(ψ, W, ft * dx, F, F̃)

    return ψ
end # T4S

"""
    T6S(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T4S`](@ref)
"""
function T6S(ψ, W, dx, F, F̃)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T4S(ψ, W, ft * dx, F, F̃)
    ψ = T4S(ψ, W, bt * dx, F, F̃)
    ψ = T4S(ψ, W, ft * dx, F, F̃)

    return ψ
end #T6S

"""
    T8S(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T6S`](@ref)
"""
function T8S(ψ, W, dx, F, F̃)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T6S(ψ, W, ft * dx, F, F̃)
    ψ = T6S(ψ, W, bt * dx, F, F̃)
    ψ = T6S(ψ, W, ft * dx, F, F̃)

    return ψ
end #T8S