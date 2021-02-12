"""
    T1A_SS!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic first order integrator of type A for the Sasa-Satsuma equation. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T1A_SS!(ψₒ, ψᵢ, dx, ops)
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

    set_u!(ops.B̂B, ψₒ)
    step!(ops.B̂B) 
    ψₒ .= ops.B̂B.u
end #T₁ʰ
