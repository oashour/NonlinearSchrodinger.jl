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

struct Sim{TT<:Real, F}
    λ::Complex{TT}
    T::TT
    Ω::TT
    box::Box{TT}
    ψ₀::Array{Complex{TT}, 1}
    T̂::F
    α::TT
    ϵ::TT
    αₚ::TT
    ψ::Array{Complex{TT}, 2}
    ψ̃::Array{Complex{TT}, 2}
    E::Array{TT, 1}
    PE::Array{TT, 1}
    KE::Array{TT, 1}
    N::Array{TT, 1}
    P::Array{TT, 1}
end # Simulation

function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}, T̂; α = 0.0, ϵ = 0.0, αₚ = 0.0) where TT <: Real
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
    sim = Sim(λ, T, Ω, box, ψ₀, T̂, α, ϵ, αₚ, ψ, ψ̃, E, PE, KE, N, P)
   return sim
end #init_sim

struct Operators{T, DispFunc, FFTPlan, InvFFTPlan}
    K̂::DispFunc
    F̂::FFTPlan
    F̃̂::InvFFTPlan
    B̂::DiffEqBase.DEIntegrator
    B̂B::DiffEqBase.DEIntegrator
    ψ₁::T
    ψ₂::T
    ψ₃::T
end # Simulation

struct Ks
    α::Float64
    ϵ::Float64
    ω::Vector{Float64}
end
  
@memoize function (K̂::Ks)(dx)
    if K̂.α == 0 && K̂.ϵ == 0 # Cubic NLSE
        @info "Caching Cubic NLS dispersoin at dx = $dx"
        ifftshift(cis.(dx*K̂.ω.^2/2))
    elseif K̂.α != 0 || K̂.ϵ != 0 # Hirota Equation
        @info "Caching Hirota/SS dispersion at dx = $dx"
        ifftshift(cis.(dx*(K̂.ω.^2/2 .- K̂.α*K̂.ω.^3)))
    end
end

function Operators(sim)
    # Generate FFT Plans to optimize performance
    @info "Generating FFT plans"
    F̂ = plan_fft!(@view sim.ψ[:, 1]) # 26 allocs
    F̃̂ = plan_ifft!(@view sim.ψ[:, 1]) # 34 allocs

    # This stuff should be user controllable
    t_algo = BS3()
    t_order = 2
    # Create the integrator for the Burger term
    if sim.α != 0 && sim.ϵ == 0
        @info "Setting up Burger operator for Hirota equation"
        D = CenteredDifference(1, t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        function MH!(du, u,p,t)
            du .= 6*sim.α*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(MH!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
        B̂B = B̂
    elseif sim.α == 0 && sim.ϵ != 0
        D = CenteredDifference(1, t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        function MS1!(du, u,p,t)
            du .= 6*sim.ϵ*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(MS1!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
        function MS2!(du, u,p,t)
            du .= 3*sim.ϵ*D*Q*abs2.(u).*u
        end
        prob = ODEProblem(MS2!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂B = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
    else
        function M0!(du, u,p,t)
            du .= 0*u 
        end
        prob = ODEProblem(M0!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false)  
    end
    
    ψ₁ = similar(sim.ψ₀)
    ψ₂ = similar(sim.ψ₀)
    ψ₃ = similar(sim.ψ₀)
    ops = Operators(Ks(sim.α, sim.ϵ, sim.box.ω), F̂, F̃̂, B̂, B̂B, ψ₁, ψ₂, ψ₃)

    return ops
end