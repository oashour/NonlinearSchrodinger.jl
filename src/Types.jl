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
    αₚ::TT
    ψ::Array{Complex{TT}, 2}
    ψ̃::Array{Complex{TT}, 2}
    E::Array{TT, 1}
    PE::Array{TT, 1}
    KE::Array{TT, 1}
    N::Array{TT, 1}
    P::Array{TT, 1}
end # Simulation

function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}, T̂; α = 0.0, αₚ = 0.0) where TT <: Real
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
    sim = Sim(λ, T, Ω, box, ψ₀, T̂, α, αₚ, ψ, ψ̃, E, PE, KE, N, P)
   return sim
end #init_sim

struct Calc{TT<:Real}
    λ::Vector{Complex{TT}}
    T::Vector{Complex{TT}}
    Ω::Vector{Complex{TT}}
    χ::Vector{Complex{TT}}
    tₛ::Vector{TT}
    xₛ::Vector{TT}
    seed::String # Should be an enum
    box::Box{TT}
    m::TT
    ψ::Matrix{Complex{TT}}
    ψ̃::Matrix{Complex{TT}}
    E::Vector{TT}
    PE::Vector{TT}
    KE::Vector{TT}
    N::Vector{TT}
    P::Vector{TT}
end # Simulation

function Calc(λ::Array{Complex{TT}}, tₛ, xₛ, seed, box; m=0.0) where TT <: Real
    if ~(length(λ) == length(tₛ) == length(xₛ))
        throw(ArgumentError("Length of shifts and eigenvalue array should be the same."))
    end
    # Elliptic modulus
    @assert m <= 1 && m >= 0
    Ω = similar(λ)
    T = similar(Ω)
    χ = similar(Ω)
    if seed == "0"
        Ω = zeros(Complex{Float64}, length(λ))
        @. χ = 0.5*acos(Ω/2)
    elseif seed == "exp"
        @. Ω = 2*sqrt(1 + λ^2)
        @. χ = 0.5*acos(Ω/2)
    elseif seed == "dn"
        @. Ω = 2*sqrt(1 + (λ - m/4/λ)^2)
        @. χ = 0.5*acos(Ω/2)
    elseif seed == "cn"
        @. Ω = 2*sqrt(m)*sqrt(1 + 1/m*(λ - 1/4/λ)^2)
        @. χ = 0.5*acos(Ω/2/sqrt(m))
    end
    @. T = 2π/Ω
    #elseif strcmp(seed, 'dn(t;k)')
    #    kappa = sqrt(1+(l - k^2/4./l).^2);  % Half the principal wave number
    #    chi = 0.5*acos(kappa);
    #elseif strcmp(seed, 'cn(t;k)')
    #    kappa = k*sqrt(1+1/k^2*(l - 1/4./l).^2); 
    #    chi  = 0.5*acos(kappa/k); 
    #end

    ψ = Array{Complex{TT}}(undef, box.Nₜ, box.Nₓ)
    ψ̃ = similar(ψ)
    E = zeros(box.Nₓ)
    PE = similar(E)    
    KE = similar(E)
    N = similar(E)
    P = similar(E)

    calc = Calc(λ, T, Ω, χ, tₛ, xₛ, seed, box, m, ψ, ψ̃, E, PE, KE, N, P)
end #init_sim

struct Operators{T, DispFunc, FFTPlan, InvFFTPlan}
    K̂::DispFunc
    F̂::FFTPlan
    F̃̂::InvFFTPlan
    B̂::DiffEqBase.DEIntegrator
    ψ₁::T
    ψ₂::T
    ψ₃::T
end # Simulation

struct Ks
    α :: Float64
    ω :: Vector{Float64}
end

@memoize function K_cubic(dx::Real, ω)
    ifftshift(cis.(dx*ω.^2/2))
end
@memoize function K_hirota(dx::Real, ω, α)
    ifftshift(cis.(dx*(ω.^2/2 .- α*ω.^3)))
end

function (K̂::Ks)(dx)
    if K̂.α == 0
        #ifftshift(cis.(dx*K̂.ω.^2/2))
        K_cubic(dx, K̂.ω)
    else
        #ifftshift(cis.(dx*(K̂.ω.^2/2 - K̂.α*K̂.ω.^3)))
        K_hirota(dx, K̂.ω, K̂.α)
    end
end

function Operators(sim)
    # Generate FFT Plans to optimize performance
    @info "Generating FFT plans"
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

    # This stuff should be user controllable
    t_algo = BS3()
    t_order = 2
    # Create the integrator for the Burger term
    if sim.α > 0
        D = CenteredDifference(1, t_order, sim.box.dt, sim.box.Nₜ) 
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
    
    ψ₁ = similar(sim.ψ₀)
    ψ₂ = similar(sim.ψ₀)
    ψ₃ = similar(sim.ψ₀)
    ops = Operators(Ks(sim.α, sim.box.ω), F̂, F̃̂, B̂, ψ₁, ψ₂, ψ₃)

    return ops
end