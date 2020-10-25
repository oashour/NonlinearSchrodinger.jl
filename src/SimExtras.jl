export Sim, Box
export params, print
export ψ₀_periodic

###########################################################################
# Types
###########################################################################
struct Box{TT<:Real}
    t::Array{TT, 1}
    ω::Array{TT, 1}
    x::Array{TT, 1}
    Nₜ::Int64
    Nₓ::Int64
    dt::TT
    dx::TT
    n_periods::Int64
end #SimulationBox

function Box(xᵣ::Pair, T; dx = 1e-3, Nₜ = 256, n_periods = 1)
    println("==========================================")
    println("Initializing simulation box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ.")
    T = n_periods * T
    println("Longitudinal range is [$(xᵣ.first), $(xᵣ.second)], transverse range is [$(-T/2), $(T/2))")
    dt = T / Nₜ
    t = dt * collect((-Nₜ/2:Nₜ/2-1))

    x = collect(xᵣ.first:dx:xᵣ.second)
    Nₓ = length(x)
    ω = 2π/T * collect((-Nₜ/2:Nₜ/2-1))

    box = Box(t, ω, x, Nₜ, Nₓ, dt, dx, n_periods)

    println("Done computing t, x, ω")
    println("==========================================")

    return box
end

mutable struct Sim{TT<:Real}
    λ::Complex{TT}
    T::TT
    Ω::TT
    box::Box{TT}
    ψ₀::Array{Complex{TT}, 1}
    t_order::Int64
    x_order::Int64
    α::TT
    αₚ::TT
    solved::Bool
    ψ::Array{Complex{TT}, 2}
    spectrum_computed::Bool
    ψ̃::Array{Complex{TT}, 2}
    energy_computed::Bool
    E::Array{TT, 1}
    PE::Array{TT, 1}
    KE::Array{TT, 1}
    dE::Array{TT, 1}
    N::Array{TT, 1}
    P::Array{TT, 1}
end # Simulation

function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}; t_order = 2, x_order = 2, α = 0.0, αₚ = 0.0) where TT <: Real
    ψ = Array{Complex{TT}}(undef, box.Nₜ, box.Nₓ)
    ψ̃ = similar(ψ)
    E = zeros(box.Nₓ)
    PE = similar(E)    
    KE = similar(E)
    dE = similar(E)
    N = similar(E)
    P = similar(E)
    if αₚ < 0.0 
        throw(ArgumentError("αₚ < 0. Set αₚ = 0 to disable pruning or αₚ = Inf for fixed pruning. See documentation)"))
    end
    if αₚ > 0.0 && box.n_periods == 1 
        throw(ArgumentError("Pruning is only applicable when n_periods > 1."))
    end
    # Compute some parameters
    λ, T, Ω = params(λ = λ)
    sim = Sim(λ, T, Ω, box, ψ₀, t_order, x_order, α, αₚ, false, ψ, false, ψ̃, false, E, PE, KE, dE, N, P)
   return sim
end #init_sim

###########################################################################
# Utility
###########################################################################
function params(; kwargs...)
    if length(kwargs) != 1
        throw(ArgumentError("You have either specified too few or too many parameters. You must specify one and only one of the following options: λ, Ω, T, a."))
    end
    param = Dict(kwargs)
    if :a in keys(param)
        λ = im * sqrt(2 * param[:a])
        T = π/sqrt(1 - imag(λ)^2)
        Ω = 2π/T
        println("Passed a = $(param[:a]), computed λ = $λ, T = $T and Ω = $Ω")
    elseif :λ  in keys(param)
        λ = param[:λ]
        T = π/sqrt(1 - imag(λ)^2)
        Ω = 2π/T
        println("Passed λ=$λ, computed T = $T and Ω = $Ω")
    elseif :Ω in keys(param)
        λ = im * sqrt((1 - (param[:Ω] / 2)^2))
        T = π/sqrt(1 - imag(λ)^2)
        Ω = param[:Ω]
        println("Passed Ω=$Ω, computed λ = $λ and T = $T")
    elseif :T in keys(param)
        λ = im * sqrt((1 - ((2*π/param[:T]) / 2)^2))
        T = param[:T]
        Ω = 2π/T
        println("Passed T = $T, computed λ = $λ and Ω = $Ω")
    end

    return λ, T, Ω
end #compute_parameters

"""
    function print(sim::Sim)

Prints information about the `Sim` instance `sim`
"""
function print(sim::Sim)
    println("Box Properties:")
    println("------------------------------------------")
    println("dx = $(sim.box.dx) (Nₓ = $(sim.box.Nₓ))")
    println("Nₜ = $(sim.box.dx) (dt = $(sim.box.dt))")
    println("Parameters:")
    println("------------------------------------------")
    println("λ = $(sim.λ)")
    println("Ω = $(sim.Ω)")
    println("T = $(sim.T)")
    println("------------------------------------------")
    # Should add information about ψ₀
    if sim.α == 0
        println("Equation: Cubic NLSE")
    elseif sim.α > 0
        println("Equation: Hirota Equation with α = $(sim.α)")
    else
        println("Unknown equation with α = $(sim.α)")
    end
    println("Algorithm of order $(sim.x_order) in x and $(sim.t_order) in t")
    println("------------------------------------------")
end

###########################################################################
# ψ₀
###########################################################################
"""
    function ψ₀_periodic(coeff, box::SimBox, params::SimParamseters; phase=0)

Computes an initial wavefunction for the `SimBox` `box`, with fundamental frequency `sim.Ω`
and coefficients ``A_1...n`` = `coeff` and an overall phase `exp(i phase t)`, i.e. of the form:

**TODO**: insert latex form here

See also: [`init_sim`](@ref)
"""
function ψ₀_periodic(coeff::Array, box::Box, Ω; phase=0)
    println("==========================================")
    println("Initializing periodic ψ₀")
    for (n, An) in enumerate(coeff)
        if abs(An) >= 1
            println("The absolute value of the coefficient A($(n+1)) = $(An) is greater than 1. psi_0 needs to be normalizable.")
        else
            println("A($(n)) = $(An)")
        end #if
    end #for

    println("Computing A₀ to preserve normalization.")
    A0 = sqrt(1 - 2 * sum([abs(An) .^ 2 for An in coeff]))
    println("Computed A₀ = $A0")
    #A0 = 1
    if phase != 0
        str = "ψ₀ = exp(i $phase t) ($A0 + "
    else
        str = "ψ₀ = $A0 + "
    end

    ψ₀ = A0 * ones(box.Nₜ)
    for (n, An) in enumerate(coeff)
        ψ₀ += 2 * An * cos.(n * Ω * box.t)
        str = string(str, "2 × $An × cos($n × $(Ω) t)")
    end #for
    # Multiply by the overall phase
    ψ₀ = exp.(im * phase * box.t) .* ψ₀

    if phase != 0
        str = string(str, ")")
    end

    println(str)
    println("==========================================")

    return ψ₀
end #psi0_periodic