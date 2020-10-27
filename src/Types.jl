@enum AlgorithmVariant A B
@enum AlgorithmType TripleJump SuzukiFractal MultiProduct Optimized 

struct Box{TT<:Real}
    t::Array{TT, 1}
    ω::Array{TT, 1}
    x::Array{TT, 1}
    Nₜ::Int64
    Nₓ::Int64
    dt::TT
    dx::TT
    n_periods::Int64
end 

function Box(xᵣ::Pair, T; dx = 1e-3, Nₜ = 256, n_periods = 1)
    @info "Initializing simulation box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ."
    T = n_periods * T
    println("Longitudinal range is [$(xᵣ.first), $(xᵣ.second)], transverse range is [$(-T/2), $(T/2))")
    dt = T / Nₜ
    t = dt * collect((-Nₜ/2:Nₜ/2-1))

    x = collect(xᵣ.first:dx:xᵣ.second)
    Nₓ = length(x)
    ω = 2π/T * collect((-Nₜ/2:Nₜ/2-1))

    box = Box(t, ω, x, Nₜ, Nₓ, dt, dx, n_periods)

    @info "Done computing t, x, ω"

    return box
end

struct Algorithm{F}
    x_order::Int64
    t_order::Int64
    variant::AlgorithmVariant
    type::AlgorithmType
    name::String
    T̂::F
end #SimulationBox

function Algorithm(;α = 0, x_order = 2, t_order = 2, variant=A, type=TripleJump)
    if order <= 2 && α == 0
        @info "Searching for algorithm of order $x_order in x, variant $variant"
    elseif order > 2 && α == 0
        @info "Searching for algorithm of order $x_order in x, variant $variant and type $type"
    elseif order <= 2 && α >= 0
        @info "Searching for algorithm of order $x_order in x, $t_order in t, variant $variant"
    elseif order > 2 && α >= 0
        @info "Searching for algorithm of order $x_order in x, $t_order in t, variant $variant"
    end

    return box
end

struct Sim{TT<:Real}
    λ::Complex{TT}
    T::TT
    Ω::TT
    box::Box{TT}
    ψ₀::Array{Complex{TT}, 1}
    t_order::Int64
    x_order::Int64
    variant::String #Should be that stuff like :A
    α::TT
    αₚ::TT
    ψ::Array{Complex{TT}, 2}
    ψ̃::Array{Complex{TT}, 2}
    E::Array{TT, 1}
    PE::Array{TT, 1}
    KE::Array{TT, 1}
    N::Array{TT, 1}
    P::Array{TT, 1}
end # Simulation

function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}; t_order = 2, x_order = 2, variant = "A", α = 0.0, αₚ = 0.0) where TT <: Real
    ψ = Array{Complex{TT}}(undef, box.Nₜ, box.Nₓ)
    ψ̃ = similar(ψ)
    E = zeros(box.Nₓ)
    PE = similar(E)    
    KE = similar(E)
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
    sim = Sim(λ, T, Ω, box, ψ₀, t_order, x_order, variant, α, αₚ, ψ, ψ̃, E, PE, KE, N, P)
   return sim
end #init_sim
struct Operators{NLSIntegrator, DispFunc, FFTPlan, InvFFTPlan}
    T̂::NLSIntegrator
    K̂::DispFunc
    F̂::FFTPlan
    F̃̂::InvFFTPlan
    B̂::DiffEqBase.DEIntegrator
end # Simulation

function Operators(sim)
    # Generate FFT Plans to optimize performance
    #println("Generating FFT Plan")
    F̂ = plan_fft!(@view sim.ψ[:, 1]) # 26 allocs
    F̃̂ = plan_ifft!(@view sim.ψ[:, 1]) # 34 allocs

    # Cache the kinetic factor
    function K̂(α)
        fun = if α == 0 
            @memoize function K_cubic(dx::Real)
                @debug "Computing and caching K(dx = $dx) for cubic NLSE"
                ifftshift(cis.(dx*sim.box.ω.^2/2))
            end
            K_cubic
        elseif α > 0
            @memoize function K_hirota(dx::Real)
                @debug "Computing and caching K(dx = $dx) for Hirota Equation"
                ifftshift(cis.(dx*(sim.box.ω.^2/2 - sim.α*sim.box.ω.^3)))
            end
            K_hirota
        end
        return fun
    end

    # Select algorithm
    t_algo = BS3() # Default Algorithm for t
    if sim.x_order == 1 && sim.α == 0 && sim.variant == "A"
        T̂ = T1A
    elseif sim.x_order == 1 && sim.α == 0 && sim.variant == "B"
        T̂ = T1B
    elseif sim.x_order == 2 && sim.α == 0 && sim.variant == "A"
        T̂ = T2A
    elseif sim.x_order == 2 && sim.α == 0 && sim.variant == "B"
        T̂ = T2B
    elseif sim.x_order == 4 && sim.α == 0 && sim.variant == "A"
        T̂ = T4A_TJ
    elseif sim.x_order == 4 && sim.α == 0 && sim.variant == "FA"
        T̂ = T4A_SF
    elseif sim.x_order == 4 && sim.α == 0 && sim.variant == "B"
        T̂ = T4B_TJ
    elseif sim.x_order == 6 && sim.α == 0 && sim.variant == "A"
        T̂ = T6A_TJ
    elseif sim.x_order == 6 && sim.α == 0 && sim.variant == "B"
        T̂ = T6B_TJ
    elseif sim.x_order == 8 && sim.α == 0 && sim.variant == "A"
        T̂ = T8A_TJ
    elseif sim.x_order == 8 && sim.α == 0 && sim.variant == "B"
        T̂ = T8B_TJ
    elseif sim.x_order == 1 && sim.α >= 0
        T̂ = T1A_H 
    elseif sim.x_order == 2 && sim.α >= 0
        T̂ = T2A_H
    elseif sim.x_order == 4 && sim.α >= 0
        T̂ = T4A_H 
        t_algo = BS3()
    else
        throw(ArgumentError("No solver available for order $x_order in x and order $t_order in t with α=$α"))
    end

    # Create the integrator for the Burger term
    if sim.α > 0
        D = CenteredDifference(1, sim.t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        function MB!(du, u,p,t)
            du .= 6*sim.α*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(MB!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
    else
        function M0!(du, u,p,t)
            du .= 0*u 
        end
        prob = ODEProblem(M0!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
    end

    ops = Operators(T̂, K̂(sim.α), F̂, F̃̂, B̂)

    return ops
end