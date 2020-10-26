"""
    T₁ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₁ˢ(ψ, dx, ops)

    # Nonlinear
    @. ψ = cis(dx * (-1*abs2(ψ)))*ψ

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx) .* ψ
    ops.F̃̂*ψ
     
end #T₁ʰ
"""
    T₂ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T₂ˢ(ψ, dx, ops)
    # Nonlinear
    @. ψ = cis(dx/2 * (-1*abs2(ψ)))*ψ

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx) .* ψ
    ops.F̃̂*ψ

    # Nonlinear
    @. ψ = cis(dx/2 * (-1*abs2(ψ)))*ψ

    return ψ
end #T₂ˢ

"""
    T₄ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T₄ˢ(ψ, dx, ops)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₂ˢ(ψ, ft*dx, ops)
    ψ = T₂ˢ(ψ, bt*dx, ops)
    ψ = T₂ˢ(ψ, ft*dx, ops)

    return ψ
end # T₄ˢ

"""
    T₆ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₄ˢ`](@ref)
"""
function T₆ˢ(ψ, dx, ops)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₄ˢ(ψ, ft*dx, ops)
    ψ = T₄ˢ(ψ, bt*dx, ops)
    ψ = T₄ˢ(ψ, ft*dx, ops)

    return ψ
end #T6S

"""
    T₈ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₆ˢ`](@ref)
"""
function T₈ˢ(ψ, dx, ops)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T₆ˢ(ψ, ft*dx, ops)
    ψ = T₆ˢ(ψ, bt*dx, ops)
    ψ = T₆ˢ(ψ, ft*dx, ops)

    return ψ
end #T8S