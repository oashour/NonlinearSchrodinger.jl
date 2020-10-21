export Algorithm

@enum Algorithm A2 A4S A6S A8S

"""
    solve!(sim::Simulation)

Solves the `Simulation` object `sim` using the techniques its attributes specify.

See also: [`init_sim`](@ref), [`NLSS.Plotter.plot_ψ`](@ref)
"""
function solve!(sim::Sim)
    #println("==========================================")
    #println("Solving cubic NLSE with the following options:")
    # Copy in x = 0 array
    sim.ψ[1, :] = sim.ψ₀
    # Check for pruning and calculate indices
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
    end
    # Generate FFT Plans to optimize performance
    #println("Generating FFT Plan")
    F = plan_fft!(sim.ψ[1, :]) # Plan
    F̃ = plan_ifft!(sim.ψ[1, :]) # Plan
    W = ifftshift(cis.(sim.box.dx * sim.box.ω .^ 2 / 2))

    # Print info about simulation
    #print(sim)

    # Start loop and pick algorithm
    #@showprogress 1 "Computing..." for i = 1:sim.box.Nₓ-1
    for i = 1:sim.box.Nₓ-1
        if sim.algorithm === A2  
            sim.ψ[i+1, :] = T2(sim.ψ[i,:], W, sim.box.dx, F, F̃)
        elseif sim.algorithm === A4S
            sim.ψ[i+1, :] = T4S(sim.ψ[i, :], sim.box.ω, sim.box.dx, F̃)
        elseif sim.algorithm == A6S
            sim.ψ[i+1, :] = T6S(sim.ψ[i, :], sim.box.ω, sim.box.dx, F)
        elseif sim.algorithm == A8S
            sim.ψ[i+1, :] = T8S(sim.ψ[i, :], sim.box.ω, sim.box.dx, F)
        #else # not needed with enums
        #    throw(ArgumentError("Algorithm type unknown, please check the documentation"))
        end
        # Pruning
        # TODO: get rid of this extra ψ and rewrite more elegantly
        if sim.αₚ > 0 
            ψ = sim.ψ[i+1, :]
            F*ψ
            @views ψ[ind_p] .*= exp.(-sim.αₚ*abs.(ψ[ind_p]))
            inv(F)*ψ
            sim.ψ[i+1, :] = ψ
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
    for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end
    # Kinetic
    F*ψ
    for i in 1:length(W)
        ψ[i] *= W[i]
    end
    F̃*ψ

    # Nonlinear
    for i in 1:length(ψ)
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
function T4S(ψ, ω, dx, F)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2(ψ, ω, ft * dx, F)
    ψ = T2(ψ, ω, bt * dx, F)
    ψ = T2(ψ, ω, ft * dx, F)

    return ψ
end # T4S

"""
    T6S(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T4S`](@ref)
"""
function T6S(ψ, ω, dx, F)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T4S(ψ, ω, ft * dx, F)
    ψ = T4S(ψ, ω, bt * dx, F)
    ψ = T4S(ψ, ω, ft * dx, F)

    return ψ
end #T6S

"""
    T8S(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T6S`](@ref)
"""
function T8S(ψ, ω, dx, F)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T6S(ψ, ω, ft * dx, F)
    ψ = T6S(ψ, ω, bt * dx, F)
    ψ = T6S(ψ, ω, ft * dx, F)

    return ψ
end #T8S