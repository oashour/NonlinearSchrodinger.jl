"""
    T₁ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₁ʰ(ψ, K, dx, F, F̃, integrator)

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx * (-1*abs2(ψ[i]))) 
    end

    # Kinetic
    KK = K(dx)
    F*ψ 
    @inbounds for i in 1:length(KK)
        ψ[i] *= KK[i]
    end
    F̃*ψ 

    # Burger
    set_u!(integrator, ψ) #1.14k allocs
    step!(integrator) #1.14k allocs
    ψ .= integrator.u
end #T₁ʰ
"""
    T₂ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₂ʰ(ψ, K, dx, F, F̃, integrator)
    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    # Kinetic
    KK = K(dx/2)
    F*ψ 
    @inbounds for i in 1:length(ψ)
        ψ[i] *= KK[i]
    end
    F̃*ψ 

    # Burger
    set_u!(integrator, ψ)
    step!(integrator)
    ψ = integrator.u

    # Kinetic
    F*ψ 
    @inbounds for i in 1:length(ψ)
        ψ[i] *= KK[i]
    end
    F̃*ψ 

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    return ψ
end #T₁ʰ