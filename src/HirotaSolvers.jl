"""
    T₁ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T1A_H!(ψₒ, ψᵢ, dx, ops)
    # Nonlinear
    @. ψₒ = cis(dx * (-1*abs2(ψᵢ)))*ψᵢ

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx) .* ψₒ
    ops.F̃̂*ψₒ

    # Burger
    set_u!(ops.B̂, ψₒ)
    step!(ops.B̂) 
    ψₒ .= ops.B̂.u
end #T₁ʰ
"""
    T₂ʰ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T2A_H!(ψₒ, ψᵢ, dx, ops)
    # Nonlinear
    @. ψₒ = cis(dx/2 * (-1*abs2(ψᵢ)))*ψᵢ

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx/2) .* ψₒ # 3 allocs
    ops.F̃̂*ψₒ

    # Burger
    set_u!(ops.B̂, ψₒ) # 0 allocs 
    step!(ops.B̂) # 3-4 allocs
    ψₒ .= ops.B̂.u # 2 allocs

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx/2) .* ψₒ # 3 allocs
    ops.F̃̂*ψₒ

    # Nonlinear
    @. ψₒ = cis(dx/2 * (-1*abs2(ψₒ)))*ψₒ
end 