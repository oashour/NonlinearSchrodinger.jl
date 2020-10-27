"""
    T₁ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T1A(ψ, dx, ops)

    # Nonlinear
    @. ψ = cis(dx * (-1*abs2(ψ)))*ψ

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx) .* ψ
    ops.F̃̂*ψ
     
end

"""
    T₁ᵇ⁽ˢ⁾(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T1B(ψ, dx, ops)

    # Dispersion
    ψ .= ops.K̂(dx) .* ψ

    # Nonlinear
    ops.F̃̂*ψ
    ψ = cis.(dx .* (-1*abs2.(ψ))).*ψ
    ops.F̂*ψ

end

"""
    function T2A(ψ, dx, ops)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T2A(ψ, dx, ops)
    # Nonlinear
    @. ψ = cis(dx/2 * (-1*abs2(ψ)))*ψ

    # Dispersion
    ops.F̂*ψ 
    ψ .= ops.K̂(dx) .* ψ
    ops.F̃̂*ψ

    # Nonlinear
    @. ψ = cis(dx/2 * (-1*abs2(ψ)))*ψ

    return ψ
end 

"""
    T2B(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic second order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref)
"""
function T2B(ψ, dx, ops)

    # Dispersion
    ψ .= ops.K̂(dx/2) .* ψ

    # Nonlinear
    ops.F̃̂*ψ
    ψ = cis.(dx .* (-1*abs2.(ψ))).*ψ
    ops.F̂*ψ

    # Dispersion
    ψ .= ops.K̂(dx/2) .* ψ

end

"""
    T4A_TJ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T4A_TJ(ψ, dx, ops)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2A(ψ, ft*dx, ops)
    ψ = T2A(ψ, bt*dx, ops)
    ψ = T2A(ψ, ft*dx, ops)

    return ψ
end

"""
    T4B_TJ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T4B_TJ(ψ, dx, ops)
    s = 2^(1 / 3)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T2B(ψ, ft*dx, ops)
    ψ = T2B(ψ, bt*dx, ops)
    ψ = T2B(ψ, ft*dx, ops)

    return ψ
end

"""
    T₆ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₄ˢ`](@ref)
"""
function T₆ᵃ⁽ˢ⁾(ψ, dx, ops)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T4A_TJ(ψ, ft*dx, ops)
    ψ = T4A_TJ(ψ, bt*dx, ops)
    ψ = T4A_TJ(ψ, ft*dx, ops)
FT4A_TJ
    return ψ
end

"""
    T₆ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₄ˢ`](@ref)
"""
function T6B_TJ(ψ, dx, ops)

    s = 2^(1 / 5)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T4B_TJ(ψ, ft*dx, ops)
    ψ = T4B_TJ(ψ, bt*dx, ops)
    ψ = T4B_TJ(ψ, ft*dx, ops)

    return ψ
end

"""
    T₈ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₆ˢ`](@ref)
"""
function T8A_TJ(ψ, dx, ops)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T6A_TJ(ψ, ft*dx, ops)
    ψ = T6A_TJ(ψ, bt*dx, ops)
    ψ = T6A_TJ(ψ, ft*dx, ops)

    return ψ
end 

"""
    T₈ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic eighth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₆ˢ`](@ref)
"""
function T8B_TJ(ψ, dx, ops)

    s = 2^(1 / 7)
    os = 1 / (2 - s)

    ft = os
    bt = -s * os

    ψ = T6B_TJ(ψ, ft*dx, ops)
    ψ = T6B_TJ(ψ, bt*dx, ops)
    ψ = T6B_TJ(ψ, ft*dx, ops)

    return ψ
end 

# Suzuki Fractal

"""
    T₄ˢ(ψ, ω, dx, F)
FT4A_TJ
Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T4A_SF(ψ, dx, ops)
    s = 4^(1 / 3)
    os = 1 / (4 - s)

    ft = os
    bt = -s * os

    ψ = T2A(ψ, ft*dx, ops)
    ψ = T2A(ψ, ft*dx, ops)
    ψ = T2A(ψ, bt*dx, ops)
    ψ = T2A(ψ, ft*dx, ops)
    ψ = T2A(ψ, ft*dx, ops)

    return ψ
end