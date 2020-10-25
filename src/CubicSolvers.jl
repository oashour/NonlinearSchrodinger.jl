"""
    T₁ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₁ˢ(ψ, K, dx, F, F̃, integrator= 0)

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
end #T₁ʰ
"""
    T₂ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₂ˢ(ψ, K, dx, F, F̃, integrator = 0)
    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end
    # Kinetic
    KK = K(dx)
    F*ψ # 0 allocs
    @inbounds for i in 1:length(ψ)
        ψ[i] *= KK[i]
    end
    F̃*ψ # 0 allocs

    # Nonlinear
    @inbounds for i in 1:length(ψ)
        ψ[i] *= cis(dx/2 * (-1*abs2(ψ[i]))) 
    end

    return ψ
end #T₂ˢ

"""
    T₄ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T₄ˢ(ψ, K, dx, F, F̃, integrator = 0)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₂ˢ(ψ, K, ft * dx, F, F̃)
    ψ = T₂ˢ(ψ, K, bt * dx, F, F̃)
    ψ = T₂ˢ(ψ, K, ft * dx, F, F̃)

    return ψ
end # T₄ˢ

"""
    T₆ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₄ˢ`](@ref)
"""
function T₆ˢ(ψ, K, dx, F, F̃, integrator = 0)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₄ˢ(ψ, K, ft * dx, F, F̃)
    ψ = T₄ˢ(ψ, K, bt * dx, F, F̃)
    ψ = T₄ˢ(ψ, K, ft * dx, F, F̃)

    return ψ
end #T6S

"""
    T₈ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₆ˢ`](@ref)
"""
function T₈ˢ(ψ, K, dx, F, F̃, integrator = 0)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₆ˢ(ψ, K, ft * dx, F, F̃)
    ψ = T₆ˢ(ψ, K, bt * dx, F, F̃)
    ψ = T₆ˢ(ψ, K, ft * dx, F, F̃)

    return ψ
end #T8S