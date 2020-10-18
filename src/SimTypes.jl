#export SimulationParameters, SimulationBox, Simulation 
export init_box, init_sim, compute_parameters
export Sim, SimBox, SimParameters

struct Box
    t::Array{Float64, 1}
    ω::Array{Float64, 1}
    x::Array{Float64, 1}
    Nₜ::Int64
    Nₓ::Int64
    dt::Float64
    dx::Float64
    n_periods::Int
end #SimulationBox

mutable struct Sim
    λ::Complex{Float64}
    T::Float64
    Ω::Float64
    box::Box
    ψ₀::Array{Complex{Float64}, 1}
    algorithm::String
    αₚ::Float64
    solved::Bool
    ψ::Array{Complex{Float64}, 2}
    spectrum_computed::Bool
    ψ̃::Array{Complex{Float64}, 2}
    energy_computed::Bool
    E::Array{Float64, 1}
    PE::Array{Float64, 1}
    KE::Array{Float64, 1}
    dE::Array{Float64, 1}
    N::Array{Float64, 1}
    P::Array{Float64, 1}
end # Simulation

function compute_parameters(; kwargs...)
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

function init_box(xᵣ::Pair, T::Float64; dx = 1e-3, Nₜ = 256, n_periods = 1)
    println("==========================================")
    println("Initializing simulation box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ.")
    T = n_periods * T
    println("Longitudinal range is [$(xᵣ.first), $(xᵣ.second)], transverse range is [$(-T/2), $(T/2))")
    dt = T / Nₜ
    t = dt * (-Nₜ/2:Nₜ/2-1)

    x = xᵣ.first:dx:xᵣ.second
    Nₓ = length(x)
    ω = 2π/T * (-Nₜ/2:Nₜ/2-1)

    box = Box(t, ω, x, Nₜ, Nₓ, dt, dx, n_periods)

    println("Done computing t, x, ω")
    println("==========================================")

    return box
end

function init_sim(λ, box::Box, ψ₀::Array; algorithm = "2S", αₚ = 0)
    ψ = Array{Complex{Float64}}(undef, box.Nₓ, box.Nₜ)
    ψ̃ = similar(ψ)
    E = zeros(box.Nₓ)
    PE = similar(E)    
    KE = similar(E)
    dE = similar(E)
    N = similar(E)
    P = similar(E)
    if αₚ < 0 
        throw(ArgumentError("αₚ < 0. Set αₚ = 0 to disable pruning or αₚ = Inf for fixed pruning. See documentation)"))
    end
    if αₚ > 0 && box.n_periods == 1 
        throw(ArgumentError("Pruning is only applicable when n_periods > 1."))
    end
    # Compute some parameters
    λ, T, Ω = compute_parameters(λ = λ)
    sim = Sim(λ, T, Ω, box, ψ₀, algorithm, αₚ, false, ψ, false, ψ̃, false, E, PE, KE, dE, 
    N, P)
   return sim
end #init_sim