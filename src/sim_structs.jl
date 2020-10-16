#export SimulationParameters, SimulationBox, Simulation 
export init_sim_box, init_sim_params, init_sim

struct SimulationBox
    t::Array{Float64, 1}
    ω::Array{Float64, 1}
    x::Array{Float64, 1}
    Nₜ::Int64
    Nₓ::Int64
    dt::Float64
    dx::Float64
    n_periods::Int
end #SimulationBox

struct SimulationParameters
    λ::Complex{Float64}
    Ω::Float64
    a::Float64
    T::Float64
end #SimulationParameters

mutable struct Simulation
    box::SimulationBox
    params::SimulationParameters
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

function init_sim_params(; kwargs...)
    println("==========================================")
    println("Initializing simulation parameters.")
    if length(kwargs) != 1
        throw(ArgumentError("You have either specified too few or too many parameters. You must specify one and only one of the following options: Ω, T, a, λ."))
    end
    param = Dict(kwargs)
    if :a in keys(param)
        a = param[:a]
        println("Passed auxillary eigenvalue a = $a")
        λ = im * sqrt(2 * a)
        Ω = 2 * sqrt(1 - 2 * a)
        T = 2π / Ω
    elseif :λ in keys(param)
        λ = param[:λ]
        println("Passed eigenvalue λ=$λ")
        a = imag(λ)^2 / 2
        Ω = 2 * sqrt(1 - 2 * a)
        T = 2π / Ω
    elseif :Ω in keys(param)
        Ω = param[:Ω]
        println("Passed frequency Ω=$Ω")
        a = (1 - (Ω / 2)^2) / 2
        λ = im * sqrt(2 * a)
        T = 2π / Ω
    elseif :T in keys(param)
        T = param[:T]
        println("Passed period T = $T")
        Ω = 2π / T
        a = (1 - (Ω / 2)^2) / 2
        λ = im * sqrt(2 * a)
    end
    params = SimulationParameters(λ, Ω, a, T)

    println("Calculated Result: a = $a, λ = $λ, Ω = $Ω, T =$T")
    println("==========================================")

    return params
end #init_sim_params

function init_sim_box(xᵣ::Pair, params::SimulationParameters; dx = 1e-3, Nₜ = 256, n_periods::Int = 1)
    println("==========================================")
    println("Initializing simulation box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ.")
    println("Box range is [xᵣ.first, xᵣ.second]")
    T = n_periods * params.T
    dt = T / Nₜ
    t = dt * (-Nₜ/2:Nₜ/2-1)

    x = xᵣ.first:dx:xᵣ.second
    Nₓ = length(x)
    ω = 2 * π / T * (-Nₜ/2:Nₜ/2-1)

    box = SimulationBox(t, ω, x, Nₜ, Nₓ, dt, dx, n_periods)

    println("Done computing t, x, ω")
    println("==========================================")

    return box
end

function init_sim(box::SimulationBox, params::SimulationParameters, ψ₀::Array; algorithm = "2S", αₚ = 0)
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
    #set_num_threads(1)
    sim = Simulation(box, params, ψ₀, algorithm, αₚ, false, ψ, false, ψ̃, false, E, PE, KE, dE, N, P)
   return sim
end #init_sim