"""
    T₁ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T1A_H(ψ, dx, ops)

    # Nonlinear
    @. ψ = cis(dx * (-1*abs2(ψ)))*ψ

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx) .* ψ
    ops.F̃̂*ψ

    # Burger
    set_u!(ops.B̂, ψ) #1.14k allocs
    step!(ops.B̂) #1.14k allocs
    ψ .= ops.B̂.u
end #T₁ʰ
"""
    T₂ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T2A_H(ψ, dx, ops)
    # Nonlinear
    @. ψ = cis(dx/2 * (-1*abs2(ψ)))*ψ

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx/2) .* ψ # 3 allocs
    ops.F̃̂*ψ

    # Burger
    set_u!(ops.B̂, ψ) # 0 allocs 
    step!(ops.B̂) # 3-4 allocs
    ψ .= ops.B̂.u # 2 allocs

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx/2) .* ψ # 3 allocs
    ops.F̃̂*ψ

    # Nonlinear
    @. ψ = cis(dx/2 * (-1*abs2(ψ)))*ψ

    return ψ
end #T₂ʰ
"""
    T₄ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T4A_TJ_H(ψ, dx, ops)
    # This algorithm is broken at the moment
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2A_H(ψ, ft*dx, ops)
    ψ = T2A_H(ψ, bt*dx, ops)
    ψ = T2A_H(ψ, ft*dx, ops)

    return ψ
end # T₄ˢ