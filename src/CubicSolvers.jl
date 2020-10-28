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
    @. ψ = cis(dx * (-1*abs2(ψ)))*ψ
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
    @. ψ = cis(dx * (-1*abs2(ψ)))*ψ
    ops.F̂*ψ

    # Dispersion
    ψ .= ops.K̂(dx/2) .* ψ

end
####################################################################
# Triple Jump
####################################################################
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

    ψ .= T2A(T2A(T2A(ψ, ft*dx, ops),bt*dx,ops),ft*dx,ops)

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

    ψ .= T2B(T2B(T2B(ψ, ft*dx, ops),bt*dx,ops),ft*dx,ops)

    return ψ
end

"""
    T₆ˢ(ψ, ω, dx, F)

Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic sixth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T₄ˢ`](@ref)
"""
function T6A_TJ(ψ, dx, ops)
    s = 2^(1 / 5)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    ψ .= T4A_TJ(T4A_TJ(T4A_TJ(ψ, ft*dx, ops),bt*dx,ops),ft*dx,ops)
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

    ψ .= T4B_TJ(T4B_TJ(T4B_TJ(ψ, ft*dx, ops),bt*dx,ops),ft*dx,ops)

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

    ψ .= T6A_TJ(T6A_TJ(T6A_TJ(ψ, ft*dx, ops),bt*dx,ops),ft*dx,ops)

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

    ψ .= T6B_TJ(T6B_TJ(T6B_TJ(ψ, ft*dx, ops),bt*dx,ops),ft*dx,ops)

    return ψ
end 

####################################################################
# Suzuki Fractal
####################################################################
"""
    T₄ˢ(ψ, ω, dx, F)
FT4A_TJ
Compute `ψ'`, i.e. `ψ` advanced a step `dx` forward using a symplectic fourth order
integrator. `ψ'` is defined on an FFT grid with frequencies `ω` using an FFT plan
`F`. Do not call this explicitly and use `solve!` instead.

See also: [`solve!`](@ref), [`T2`](@ref)
"""
function T4A_SF(ψ, dx, ops)
    s = 4^(1/3)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    ψ .= T2A(T2A(T2A(T2A(T2A(ψ,ft*dx, ops),ft*dx,ops),bt*dx,ops),ft*dx,ops)
               ,ft*dx,ops)

    return ψ
end

function T4B_SF(ψ, dx, ops)
    s = 4^(1/3)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    ψ .= T2B(T2B(T2B(T2B(T2B(ψ,ft*dx, ops),ft*dx,ops),bt*dx,ops),ft*dx,ops)
               ,ft*dx,ops)

    return ψ
end

function T6A_SF(ψ, dx, ops)
    s = 4^(1/5)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    ψ .= T4A_SF(T4A_SF(T4A_SF(T4A_SF(T4A_SF(ψ,ft*dx, ops),ft*dx,ops),bt*dx,ops),ft*dx,ops)
               ,ft*dx,ops)

    return ψ
end

function T6B_SF(ψ, dx, ops)
    s = 4^(1/5)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    ψ .= T4B_SF(T4B_SF(T4B_SF(T4B_SF(T4B_SF(ψ,ft*dx, ops),ft*dx,ops),bt*dx,ops),ft*dx,ops)
               ,ft*dx,ops)

    return ψ
end

function T8A_SF(ψ, dx, ops)
    s = 4^(1/7)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    ψ .= T6A_SF(T6A_SF(T6A_SF(T6A_SF(T6A_SF(ψ,ft*dx, ops),ft*dx,ops),bt*dx,ops),ft*dx,ops)
               ,ft*dx,ops)
    return ψ
end

function T8B_SF(ψ, dx, ops)
    s = 4^(1/7)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    ψ .= T6B_SF(T6B_SF(T6B_SF(T6B_SF(T6B_SF(ψ,ft*dx, ops),ft*dx,ops),bt*dx,ops),ft*dx,ops)
               ,ft*dx,ops)

    return ψ
end
####################################################################
# Multi-Product Nystrom
####################################################################
function T4A_N(ψ, dx, ops)
    ψ .= 4/3*T2A(copy(ψ), dx,ops) - 1/3*T2A(T2A(copy(ψ),dx/2,ops),dx/2,ops)
end

function T4B_N(ψ, dx, ops)
    ψ .= 4/3*T2B(copy(ψ), dx,ops) - 1/3*T2B(T2B(copy(ψ),dx/2,ops),dx/2,ops)
end

function T6A_N(ψ, dx, ops)
    ψ .= 81/40*T2A(T2A(T2A(copy(ψ),dx/3.0,ops),dx/3.0,ops),dx/3.0,ops) +
         -16/15*T2A(T2A(copy(ψ),dx/2.0,ops),dx/2.0,ops) + 
         1/24*T2A(copy(ψ),dx,ops)
end

function T6B_N(ψ, dx, ops)
    ψ .= 81/40*T2B(T2B(T2B(copy(ψ),dx/3.0,ops),dx/3.0,ops),dx/3.0,ops) +
         -16/15*T2B(T2B(copy(ψ),dx/2.0,ops),dx/2.0,ops) + 
         1/24*T2B(copy(ψ),dx,ops)
end

function T8A_N(ψ, dx, ops)
    ψ .= 1024/315*T2A(T2A(T2A(T2A(copy(ψ),dx/4.0,ops),dx/4.0,ops),dx/4.0,ops),dx/4.0,ops) +
         -729/280*T2A(T2A(T2A(copy(ψ),dx/3.0,ops),dx/3.0,ops),dx/3.0,ops) +
         16/45*T2A(T2A(copy(ψ),dx/2.0,ops),dx/2.0,ops) + 
         -1/360*T2A(copy(ψ),dx,ops)
end

function T8B_N(ψ, dx, ops)
    ψ .= 1024/315*T2B(T2B(T2B(T2B(copy(ψ),dx/4.0,ops),dx/4.0,ops),dx/4.0,ops),dx/4.0,ops) +
         -729/280*T2B(T2B(T2B(copy(ψ),dx/3.0,ops),dx/3.0,ops),dx/3.0,ops) +
         16/45*T2B(T2B(copy(ψ),dx/2.0,ops),dx/2.0,ops) + 
         -1/360*T2B(copy(ψ),dx,ops)
end

####################################################################
# Optimized
####################################################################
function T6A_OP(ψ, dx, ops)
    γ₁ = 0.392256805238773 
    γ₂ = γ₁
    γ₃ = 0.1177866066796810
    γ₄ = γ₃
    γ₅ = -0.5888399920894384
    γ₆ = γ₅
    γ₇ = 0.6575931603419684 
    γ₈ = γ₇
    γ₉ = γ₅
    γ₁₀ = γ₅
    γ₁₁ = γ₃
    γ₁₂ = γ₃
    γ₁₃ = γ₁
    γ₁₄ = γ₁
    ψ .= T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(ψ,γ₁*dx,ops),γ₂*dx,ops),
    γ₃*dx,ops),γ₄*dx,ops),γ₅*dx,ops),γ₆*dx,ops),γ₇*dx,ops),γ₈*dx,ops),γ₉*dx,ops),γ₁₀*dx,ops)
    ,γ₁₁*dx,ops),γ₁₂*dx,ops),γ₁₃*dx,ops),γ₁₄*dx,ops) 
end

function T6B_OP(ψ, dx, ops)
    γ₁ = 0.392256805238773 
    γ₂ = γ₁
    γ₃ = 0.1177866066796810
    γ₄ = γ₃
    γ₅ = -0.5888399920894384
    γ₆ = γ₅
    γ₇ = 0.6575931603419684 
    γ₈ = γ₇
    γ₉ = γ₅
    γ₁₀ = γ₅
    γ₁₁ = γ₃
    γ₁₂ = γ₃
    γ₁₃ = γ₁
    γ₁₄ = γ₁
    ψ .= T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(ψ,γ₁*dx,ops),γ₂*dx,ops),
    γ₃*dx,ops),γ₄*dx,ops),γ₅*dx,ops),γ₆*dx,ops),γ₇*dx,ops),γ₈*dx,ops),γ₉*dx,ops),γ₁₀*dx,ops)
    ,γ₁₁*dx,ops),γ₁₂*dx,ops),γ₁₃*dx,ops),γ₁₄*dx,ops) 
end

function T8A_OP(ψ, dx, ops)
    γ₁ = 0.7416703643506129534482278
    γ₂ = -0.409100825800031593997300
    γ₃ = 0.1907547102962383799538763
    γ₄ = -0.5738624711160822666563877
    γ₅ = 0.2990641813036559238444635
    γ₆ = 0.3346249182452981837849580
    γ₇ = 0.3152930923967665966320567
    γ₈ = -0.7968879393529163540197888
    γ₉ = γ₇
    γ₁₀ = γ₆
    γ₁₁ = γ₅
    γ₁₂ = γ₄
    γ₁₃ = γ₃
    γ₁₄ = γ₂
    γ₁₅ = γ₁
    ψ .= T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(T2A(ψ,γ₁*dx,ops),γ₂*dx,ops),γ₃*dx,ops),γ₄*dx,ops),γ₅*dx,ops),γ₆*dx,ops),γ₇*dx,ops),γ₈*dx,ops),γ₉*dx,ops),γ₁₀*dx,ops),γ₁₁*dx,ops),γ₁₂*dx,ops),γ₁₃*dx,ops),γ₁₄*dx,ops),γ₁₅*dx,ops)
end

function T8B_OP(ψ, dx, ops)
    γ₁ = 0.7416703643506129534482278
    γ₂ = -0.409100825800031593997300
    γ₃ = 0.1907547102962383799538763
    γ₄ = -0.5738624711160822666563877
    γ₅ = 0.2990641813036559238444635
    γ₆ = 0.3346249182452981837849580
    γ₇ = 0.3152930923967665966320567
    γ₈ = -0.7968879393529163540197888
    γ₉ = γ₇
    γ₁₀ = γ₆
    γ₁₁ = γ₅
    γ₁₂ = γ₄
    γ₁₃ = γ₃
    γ₁₄ = γ₂
    γ₁₅ = γ₁
    ψ .= T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(T2B(ψ,γ₁*dx,ops),γ₂*dx,ops),γ₃*dx,ops),γ₄*dx,ops),γ₅*dx,ops),γ₆*dx,ops),γ₇*dx,ops),γ₈*dx,ops),γ₉*dx,ops),γ₁₀*dx,ops),γ₁₁*dx,ops),γ₁₂*dx,ops),γ₁₃*dx,ops),γ₁₄*dx,ops),γ₁₅*dx,ops)
end