export solve!

function solve!(sim::Simulation)
    """
        solve!(sim::Simulation)

        Solves the `Simulation` object `sim` using the techniques its attributes specify.

        See also: [`init_sim`](@ref), [`Simulation`](@ref)
    """
    println("==========================================")
    println("Solving cubic NLSE with the following options:")
    # Copy in x = 0 array
    sim.ψ[1, :] = sim.ψ₀
    # Check for pruning and calculate indices
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
    end
    # Generate FFT Plans to optimize performance
    println("Generating FFT Plan")
    F = plan_fft!(sim.ψ[1, :]) # Plan
    F̃ = inv(F) # Inverse plan, now cached

    # Print info about simulation
    print(sim)

    # Start loop and pick algorithm
    @showprogress 1 "Computing..." for i = 1:sim.box.Nₓ-1
        if sim.algorithm == "2S"
            sim.ψ[i+1, :] = T2(sim.ψ[i, :], sim.box.ω, sim.box.dx, F)
        elseif sim.algorithm == "4S"
            sim.ψ[i+1, :] = T4S(sim.ψ[i, :], sim.box.ω, sim.box.dx, F)
        elseif sim.algorithm == "6S"
            sim.ψ[i+1, :] = T6S(sim.ψ[i, :], sim.box.ω, sim.box.dx, F)
        elseif sim.algorithm == "8S"
            sim.ψ[i+1, :] = T8S(sim.ψ[i, :], sim.box.ω, sim.box.dx, F)
        else
            throw(ArgumentError("Algorithm type unknown, please check the documentation"))
        end
        # Pruning
        # TODO: get rid of this extra ψ and rewrite more elegantly
        if sim.αₚ > 0 
            ψ = sim.ψ[i+1, :]
            F*ψ
            @views ψ[ind_p] .*= exp.(-sim.αₚ*abs.(ψ[i+1, ind_p]))
            inv(F)*ψ
            sim.ψ[i+1, :] = ψ
        end
    end #for
    sim.solved = true

    println("Computation Done!")
    println("==========================================")

    return nothing
end #solve

function T2(ψ, ω, dx, F)
    """
        T2(ψ, ω, dx, F)

        Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
        integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
        `F`. Uses second order integrator. Private function that should not be called
        explicitly. Please use `solve`.

        See also: [`solve`](@ref)
    """
    # Nonlinear
    V = -1*abs.(ψ).^2                      
    ψ .*= exp.(-im * dx/2 * (-1*abs.(ψ).^2)) 

    # Kinetic
    F*ψ
    ψ .*= ifftshift(exp.(-im * dx * ω .^ 2 / 2)) 
    inv(F)*ψ

    # Nonlinear
    ψ .*= exp.(-im * dx/2 * (-1*abs.(ψ).^2)) 

    return ψ
end #T2

function T4S(ψ, ω, dx, F)
    """
        T4S(ψ, ω, dx, F)

        Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
        integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
        `F`. Uses second order integrator. Private function that should not be called
        explicitly. Please use `solve`.

        See also: [`solve`](@ref), [`T2`](@ref)
    """
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2(ψ, ω, ft * dx, F)
    ψ = T2(ψ, ω, bt * dx, F)
    ψ = T2(ψ, ω, ft * dx, F)

    return ψ
end # T4S

function T6S(ψ, ω, dx, F)
    """
        T6S(ψ, ω, dx, F)

        Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
        integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
        `F`. Uses second order integrator. Private function that should not be called
        explicitly. Please use `solve`.

        See also: [`solve`](@ref), [`T4S`](@ref)
    """

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T4S(ψ, ω, ft * dx, F)
    ψ = T4S(ψ, ω, bt * dx, F)
    ψ = T4S(ψ, ω, ft * dx, F)

    return ψ
end #T6S

function T8S(ψ, ω, dx, F)
    """
        T8S(ψ, ω, dx, F)

        Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
        integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
        `F`. Uses second order integrator. Private function that should not be called
        explicitly. Please use `solve`.

        See also: [`solve`](@ref), [`T6S`](@ref)
    """

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T6S(ψ, ω, ft * dx, F)
    ψ = T6S(ψ, ω, bt * dx, F)
    ψ = T6S(ψ, ω, ft * dx, F)

    return ψ
end #T8S