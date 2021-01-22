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
