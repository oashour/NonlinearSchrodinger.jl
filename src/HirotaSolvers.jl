"""
    T1A_H!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic first order integrator of type A for the Hirota equation. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T1A_H!(ψₒ, ψᵢ, dx, ops)
    # Nonlinear
    @. ψₒ = cis(-dx * (-1*abs2(ψᵢ)))*ψᵢ

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
    T2A_H!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic second order integrator of type A for the Hirota equation. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T2A_H!(ψₒ, ψᵢ, dx, ops)
    # Nonlinear
    @. ψₒ = cis(-dx/2 * (-1*abs2(ψᵢ)))*ψᵢ

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx/2) .* ψₒ 
    ops.F̃̂*ψₒ

    # Burger
    set_u!(ops.B̂, ψₒ) 
    step!(ops.B̂) 
    ψₒ .= ops.B̂.u 

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx/2) .* ψₒ 
    ops.F̃̂*ψₒ

    # Nonlinear
    @. ψₒ = cis(-dx/2 * (-1*abs2(ψₒ)))*ψₒ
end 