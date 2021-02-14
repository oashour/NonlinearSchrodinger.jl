####################################################################
# Symplectic
####################################################################

"""
    T1A!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic first order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T1A!(ψₒ, ψᵢ, dx, ops)

    # Nonlinear
    @. ψₒ = cis(-dx * (-1*abs2(ψᵢ)))*ψᵢ

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx) .* ψₒ
    ops.F̃̂*ψₒ
end

"""
    T1B!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic first order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T1B!(ψₒ, ψᵢ, dx, ops)

    # Dispersion
    ψₒ .= ops.K̂(dx) .* ψᵢ

    # Nonlinear
    ops.F̃̂*ψₒ
    @. ψₒ = cis(-dx * (-1*abs2(ψₒ)))*ψₒ
    ops.F̂*ψₒ
end

"""
    T2A!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic second order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T2A!(ψₒ, ψᵢ, dx, ops)
    # Nonlinear
    @. ψₒ = cis(-dx/2 * (-1*abs2(ψᵢ)))*ψᵢ

    # Dispersion
    ops.F̂*ψₒ
    ψₒ .= ops.K̂(dx) .* ψₒ
    ops.F̃̂*ψₒ

    # Nonlinear
    @. ψₒ = cis(-dx/2 * (-1*abs2(ψₒ)))*ψₒ
end 

"""
    T2B!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic second order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T2B!(ψₒ, ψᵢ, dx, ops)

    # Dispersion
    ψₒ .= ops.K̂(dx/2) .* ψᵢ

    # Nonlinear
    ops.F̃̂*ψₒ
    @. ψₒ = cis(-dx * (-1*abs2(ψₒ)))*ψₒ
    ops.F̂*ψₒ

    # Dispersion
    ψₒ .= ops.K̂(dx/2) .* ψₒ
end

####################################################################
# Triple Jump
####################################################################
"""
    T4A_TJ!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Triple Jump Fourth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T4A_TJ!(ψₒ, ψᵢ, dx, ops)
    s = 2^(1 / 3)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    T2A!(ψₒ, ψᵢ, ft*dx, ops)
    T2A!(ψₒ, ψₒ, bt*dx, ops)
    T2A!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T4B_TJ!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Triple Jump Fourth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T4B_TJ!(ψₒ, ψᵢ, dx, ops)
    s = 2^(1 / 3)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    T2B!(ψₒ, ψᵢ, ft*dx, ops)
    T2B!(ψₒ, ψₒ, bt*dx, ops)
    T2B!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T6A_TJ!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Triple Jump Sixth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6A_TJ!(ψₒ, ψᵢ, dx, ops)
    s = 2^(1 / 5)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    T4A_TJ!(ψₒ, ψᵢ, ft*dx, ops)
    T4A_TJ!(ψₒ, ψₒ, bt*dx, ops)
    T4A_TJ!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T6B_TJ!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Triple Jump Sixth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6B_TJ!(ψₒ, ψᵢ, dx, ops)
    s = 2^(1 / 5)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    T4B_TJ!(ψₒ, ψᵢ, ft*dx, ops)
    T4B_TJ!(ψₒ, ψₒ, bt*dx, ops)
    T4B_TJ!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T8A_TJ!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Triple Jump Eighth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8A_TJ!(ψₒ, ψᵢ, dx, ops)
    s = 2^(1 / 7)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    T6A_TJ!(ψₒ, ψᵢ, ft*dx, ops)
    T6A_TJ!(ψₒ, ψₒ, bt*dx, ops)
    T6A_TJ!(ψₒ, ψₒ, ft*dx, ops)
end 

"""
    T8A_TJ!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Triple Jump Eighth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8B_TJ!(ψₒ, ψᵢ, dx, ops)
    s = 2^(1 / 7)
    os = 1 / (2 - s)
    ft = os
    bt = -s * os

    T6B_TJ!(ψₒ, ψᵢ, ft*dx, ops)
    T6B_TJ!(ψₒ, ψₒ, bt*dx, ops)
    T6B_TJ!(ψₒ, ψₒ, ft*dx, ops)
end 

####################################################################
# Suzuki Fractal
####################################################################

"""
    T4A_SF!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Suzuki's Fractal Fourth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T4A_SF!(ψₒ, ψᵢ, dx, ops)
    s = 4^(1/3)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    T2A!(ψₒ, ψᵢ, ft*dx, ops)
    T2A!(ψₒ, ψₒ, ft*dx, ops)
    T2A!(ψₒ, ψₒ, bt*dx, ops)
    T2A!(ψₒ, ψₒ, ft*dx, ops)
    T2A!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T4B_SF!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Suzuki's Fractal Fourth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T4B_SF!(ψₒ, ψᵢ, dx, ops)
    s = 4^(1/3)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    T2B!(ψₒ, ψᵢ, ft*dx, ops)
    T2B!(ψₒ, ψₒ, ft*dx, ops)
    T2B!(ψₒ, ψₒ, bt*dx, ops)
    T2B!(ψₒ, ψₒ, ft*dx, ops)
    T2B!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T6A_SF!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Suzuki's Fractal Sixth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6A_SF!(ψₒ, ψᵢ, dx, ops)
    s = 4^(1/5)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    T4A_SF!(ψₒ, ψᵢ, ft*dx, ops)
    T4A_SF!(ψₒ, ψₒ, ft*dx, ops)
    T4A_SF!(ψₒ, ψₒ, bt*dx, ops)
    T4A_SF!(ψₒ, ψₒ, ft*dx, ops)
    T4A_SF!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T6B_SF!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Suzuki's Fractal Sixth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6B_SF!(ψₒ, ψᵢ, dx, ops)
    s = 4^(1/5)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    T4B_SF!(ψₒ, ψᵢ, ft*dx, ops)
    T4B_SF!(ψₒ, ψₒ, ft*dx, ops)
    T4B_SF!(ψₒ, ψₒ, bt*dx, ops)
    T4B_SF!(ψₒ, ψₒ, ft*dx, ops)
    T4B_SF!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T8A_SF!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Suzuki's Fractal Eighth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8A_SF!(ψₒ, ψᵢ, dx, ops)
    s = 4^(1/7)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    T6A_SF!(ψₒ, ψᵢ, ft*dx, ops)
    T6A_SF!(ψₒ, ψₒ, ft*dx, ops)
    T6A_SF!(ψₒ, ψₒ, bt*dx, ops)
    T6A_SF!(ψₒ, ψₒ, ft*dx, ops)
    T6A_SF!(ψₒ, ψₒ, ft*dx, ops)
end

"""
    T8B_SF!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a symplectic Suzuki's Fractal Eighth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8B_SF!(ψₒ, ψᵢ, dx, ops)
    s = 4^(1/7)
    os = 1/(4 - s)
    ft = os
    bt = -s*os

    T6B_SF!(ψₒ, ψᵢ, ft*dx, ops)
    T6B_SF!(ψₒ, ψₒ, ft*dx, ops)
    T6B_SF!(ψₒ, ψₒ, bt*dx, ops)
    T6B_SF!(ψₒ, ψₒ, ft*dx, ops)
    T6B_SF!(ψₒ, ψₒ, ft*dx, ops)
end

####################################################################
# Multi-Product Nystrom
####################################################################
"""
    T4A_CMP!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a Chin Multi-Product Fourth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T4A_CMP!(ψₒ, ψᵢ, dx, ops)
    T2A!(ψₒ, ψᵢ, dx/2, ops)
    T2A!(ψₒ, ψₒ, dx/2, ops)
    ψₒ .*= -1/3

    T2A!(ops.ψ₁, ψᵢ, dx,ops)
    ops.ψ₁ .*= 4/3

    ψₒ .+= ops.ψ₁
end

"""
    T4B_CMP!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a Chin Multi-Product Fourth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T4B_CMP!(ψₒ, ψᵢ, dx, ops)
    T2B!(ψₒ, ψᵢ, dx/2, ops)
    T2B!(ψₒ, ψₒ, dx/2, ops)
    ψₒ .*= -1/3

    T2B!(ops.ψ₁, ψᵢ, dx,ops)
    ops.ψ₁ .*= 4/3

    ψₒ .+= ops.ψ₁
end

"""
    T6A_CMP!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a Chin Multi-Product Sixth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6A_CMP!(ψₒ, ψᵢ, dx, ops)
    T2A!(ψₒ, ψᵢ, dx/3, ops)
    T2A!(ψₒ, ψₒ, dx/3, ops)
    T2A!(ψₒ, ψₒ, dx/3, ops)
    ψₒ .*= 81/40

    T2A!(ops.ψ₁, ψᵢ, dx/2, ops)
    T2A!(ops.ψ₁, ops.ψ₁, dx/2, ops)
    ops.ψ₁ .*= -16/15

    T2A!(ops.ψ₂, ψᵢ, dx, ops)
    ops.ψ₂ .*= 1/24

    ψₒ .+= ops.ψ₁ .+ ops.ψ₂
end

"""
    T6B_CMP!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a Chin Multi-Product Sixth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6B_CMP!(ψₒ, ψᵢ, dx, ops)
    T2B!(ψₒ, ψᵢ, dx/3, ops)
    T2B!(ψₒ, ψₒ, dx/3, ops)
    T2B!(ψₒ, ψₒ, dx/3, ops)
    ψₒ .*= 81/40

    T2B!(ops.ψ₁, ψᵢ, dx/2, ops)
    T2B!(ops.ψ₁, ops.ψ₁, dx/2, ops)
    ops.ψ₁ .*= -16/15

    T2B!(ops.ψ₂, ψᵢ, dx, ops)
    ops.ψ₂ .*= 1/24

    ψₒ .+= ops.ψ₁ .+ ops.ψ₂
end

"""
    T8A_CMP!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a Chin Multi-Product Eighth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8A_CMP!(ψₒ, ψᵢ, dx, ops)
    T2A!(ψₒ, ψᵢ, dx/4, ops)
    T2A!(ψₒ, ψₒ, dx/4, ops)
    T2A!(ψₒ, ψₒ, dx/4, ops)
    T2A!(ψₒ, ψₒ, dx/4, ops)
    ψₒ .*= 1024/315

    T2A!(ops.ψ₁, ψᵢ, dx/3, ops)
    T2A!(ops.ψ₁, ops.ψ₁, dx/3, ops)
    T2A!(ops.ψ₁, ops.ψ₁, dx/3, ops)
    ops.ψ₁ .*= -729/280

    T2A!(ops.ψ₂, ψᵢ, dx/2, ops)
    T2A!(ops.ψ₂, ops.ψ₂, dx/2, ops)
    ops.ψ₂ .*= 16/45

    T2A!(ops.ψ₃, ψᵢ, dx, ops)
    ops.ψ₃ .*= -1/360

    ψₒ .+= ops.ψ₁ .+ ops.ψ₂ .+ ops.ψ₃
end

"""
    T8B_CMP!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using a Chin Multi-Product Eighth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8B_CMP!(ψₒ, ψᵢ, dx, ops)
    T2B!(ψₒ, ψᵢ, dx/4, ops)
    T2B!(ψₒ, ψₒ, dx/4, ops)
    T2B!(ψₒ, ψₒ, dx/4, ops)
    T2B!(ψₒ, ψₒ, dx/4, ops)
    ψₒ .*= 1024/315

    T2B!(ops.ψ₁, ψᵢ, dx/3, ops)
    T2B!(ops.ψ₁, ops.ψ₁, dx/3, ops)
    T2B!(ops.ψ₁, ops.ψ₁, dx/3, ops)
    ops.ψ₁ .*= -729/280

    T2B!(ops.ψ₂, ψᵢ, dx/2, ops)
    T2B!(ops.ψ₂, ops.ψ₂, dx/2, ops)
    ops.ψ₂ .*= 16/45

    T2B!(ops.ψ₃, ψᵢ, dx, ops)
    ops.ψ₃ .*= -1/360

    ψₒ .+= ops.ψ₁ .+ ops.ψ₂ .+ ops.ψ₃
end

####################################################################
# Optimized
####################################################################
"""
    T6A_Ys7!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Yoshida's s7 Symplectic Sixth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6A_Ys7!(ψₒ, ψᵢ, dx, ops)
    γ₁ = 0.78451361047755726381949763
    γ₂ = 0.23557321335935813368479318
    γ₃ = -1.17767998417887100694641568
    γ₄ = 1.31518632068391121888424973
    γ₇ = γ₁
    γ₆ = γ₂
    γ₅ = γ₃

    T2A!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2A!(ψₒ, ψₒ, γ₂*dx, ops)
    T2A!(ψₒ, ψₒ, γ₃*dx, ops)
    T2A!(ψₒ, ψₒ, γ₄*dx, ops)
    T2A!(ψₒ, ψₒ, γ₅*dx, ops)
    T2A!(ψₒ, ψₒ, γ₆*dx, ops)
    T2A!(ψₒ, ψₒ, γ₇*dx, ops)
end

"""
    T6B_Ys7!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Yoshida's s7 Symplectic Sixth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6B_Ys7!(ψₒ, ψᵢ, dx, ops)
    γ₁ = 0.78451361047755726381949763
    γ₂ = 0.23557321335935813368479318
    γ₃ = -1.17767998417887100694641568
    γ₄ = 1.31518632068391121888424973
    γ₇ = γ₁
    γ₆ = γ₂
    γ₅ = γ₃

    T2B!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2B!(ψₒ, ψₒ, γ₂*dx, ops)
    T2B!(ψₒ, ψₒ, γ₃*dx, ops)
    T2B!(ψₒ, ψₒ, γ₄*dx, ops)
    T2B!(ψₒ, ψₒ, γ₅*dx, ops)
    T2B!(ψₒ, ψₒ, γ₆*dx, ops)
    T2B!(ψₒ, ψₒ, γ₇*dx, ops)
end

"""
    T6A_KLs9!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Kahan & Li's s9 Symplectic Sixth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6A_KLs9!(ψₒ, ψᵢ, dx, ops)
    γ₁ = 0.39216144400731413927925056
    γ₂ = 0.33259913678935943859974864
    γ₃ = -0.70624617255763935980996482
    γ₄ =  0.08221359629355080023149045
    γ₅ = 0.79854399093482996339895035
    γ₉ = γ₁
    γ₈ = γ₂
    γ₇ = γ₃
    γ₆ = γ₄

    T2A!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2A!(ψₒ, ψₒ, γ₂*dx, ops)
    T2A!(ψₒ, ψₒ, γ₃*dx, ops)
    T2A!(ψₒ, ψₒ, γ₄*dx, ops)
    T2A!(ψₒ, ψₒ, γ₅*dx, ops)
    T2A!(ψₒ, ψₒ, γ₆*dx, ops)
    T2A!(ψₒ, ψₒ, γ₇*dx, ops)
    T2A!(ψₒ, ψₒ, γ₈*dx, ops)
    T2A!(ψₒ, ψₒ, γ₉*dx, ops)
end

"""
    T6B_KLs9!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Kahan & Li's s9 Symplectic Sixth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6B_KLs9!(ψₒ, ψᵢ, dx, ops)
    γ₁ = 0.39216144400731413927925056
    γ₂ = 0.33259913678935943859974864
    γ₃ = -0.70624617255763935980996482
    γ₄ =  0.08221359629355080023149045
    γ₅ = 0.79854399093482996339895035
    γ₉ = γ₁
    γ₈ = γ₂
    γ₇ = γ₃
    γ₆ = γ₄

    T2B!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2B!(ψₒ, ψₒ, γ₂*dx, ops)
    T2B!(ψₒ, ψₒ, γ₃*dx, ops)
    T2B!(ψₒ, ψₒ, γ₄*dx, ops)
    T2B!(ψₒ, ψₒ, γ₅*dx, ops)
    T2B!(ψₒ, ψₒ, γ₆*dx, ops)
    T2B!(ψₒ, ψₒ, γ₇*dx, ops)
    T2B!(ψₒ, ψₒ, γ₈*dx, ops)
    T2B!(ψₒ, ψₒ, γ₉*dx, ops)
end

"""
    T6A_Ss14!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Suzuki's s14 Symplectic Sixth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6A_Ss14!(ψₒ, ψᵢ, dx, ops)
    γ₁ = 0.392256805238773 
    γ₂ = γ₁
    γ₁₃ = γ₁
    γ₁₄ = γ₁
    γ₃ = 0.1177866066796810
    γ₄ = γ₃
    γ₁₁ = γ₃
    γ₁₂ = γ₃
    γ₅ = -0.5888399920894384
    γ₆ = γ₅
    γ₉ = γ₅
    γ₁₀ = γ₅
    γ₇ = 0.6575931603419684 
    γ₈ = γ₇

    T2A!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2A!(ψₒ, ψₒ, γ₂*dx, ops)
    T2A!(ψₒ, ψₒ, γ₃*dx, ops)
    T2A!(ψₒ, ψₒ, γ₄*dx, ops)
    T2A!(ψₒ, ψₒ, γ₅*dx, ops)
    T2A!(ψₒ, ψₒ, γ₆*dx, ops)
    T2A!(ψₒ, ψₒ, γ₇*dx, ops)
    T2A!(ψₒ, ψₒ, γ₈*dx, ops)
    T2A!(ψₒ, ψₒ, γ₉*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₀*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₁*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₂*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₃*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₄*dx, ops)
end

"""
    T6B_Ss14!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Suzuki's s14 Symplectic Sixth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T6B_Ss14!(ψₒ, ψᵢ, dx, ops)
    γ₁ = 0.392256805238773 
    γ₂ = γ₁
    γ₁₃ = γ₁
    γ₁₄ = γ₁
    γ₃ = 0.1177866066796810
    γ₄ = γ₃
    γ₁₁ = γ₃
    γ₁₂ = γ₃
    γ₅ = -0.5888399920894384
    γ₆ = γ₅
    γ₉ = γ₅
    γ₁₀ = γ₅
    γ₇ = 0.6575931603419684 
    γ₈ = γ₇

    T2B!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2B!(ψₒ, ψₒ, γ₂*dx, ops)
    T2B!(ψₒ, ψₒ, γ₃*dx, ops)
    T2B!(ψₒ, ψₒ, γ₄*dx, ops)
    T2B!(ψₒ, ψₒ, γ₅*dx, ops)
    T2B!(ψₒ, ψₒ, γ₆*dx, ops)
    T2B!(ψₒ, ψₒ, γ₇*dx, ops)
    T2B!(ψₒ, ψₒ, γ₈*dx, ops)
    T2B!(ψₒ, ψₒ, γ₉*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₀*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₁*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₂*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₃*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₄*dx, ops)
end

"""
    T8A_Ss15!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Suzuki's s15 Symplectic Eighth order integrator of type A. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8A_Ss15!(ψₒ, ψᵢ, dx, ops)
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

    T2A!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2A!(ψₒ, ψₒ, γ₂*dx, ops)
    T2A!(ψₒ, ψₒ, γ₃*dx, ops)
    T2A!(ψₒ, ψₒ, γ₄*dx, ops)
    T2A!(ψₒ, ψₒ, γ₅*dx, ops)
    T2A!(ψₒ, ψₒ, γ₆*dx, ops)
    T2A!(ψₒ, ψₒ, γ₇*dx, ops)
    T2A!(ψₒ, ψₒ, γ₈*dx, ops)
    T2A!(ψₒ, ψₒ, γ₉*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₀*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₁*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₂*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₃*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₄*dx, ops)
    T2A!(ψₒ, ψₒ, γ₁₅*dx, ops)
end

"""
    T8B_Ss15!(ψₒ, ψᵢ, dx, ops)

Compute `ψₒ`, i.e. `ψᵢ` advanced a step `dx` forward using Suzuki's s15 Symplectic Eighth order integrator of type B. The structure `ops::Operators` contains the FFT plans and the kinetic energy operators.  

See also: [`solve!`](@ref), [`Operators`](@ref)
"""
function T8B_Ss15!(ψₒ, ψᵢ, dx, ops)
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

    T2B!(ψₒ, ψᵢ, γ₁*dx, ops)
    T2B!(ψₒ, ψₒ, γ₂*dx, ops)
    T2B!(ψₒ, ψₒ, γ₃*dx, ops)
    T2B!(ψₒ, ψₒ, γ₄*dx, ops)
    T2B!(ψₒ, ψₒ, γ₅*dx, ops)
    T2B!(ψₒ, ψₒ, γ₆*dx, ops)
    T2B!(ψₒ, ψₒ, γ₇*dx, ops)
    T2B!(ψₒ, ψₒ, γ₈*dx, ops)
    T2B!(ψₒ, ψₒ, γ₉*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₀*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₁*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₂*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₃*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₄*dx, ops)
    T2B!(ψₒ, ψₒ, γ₁₅*dx, ops)
end