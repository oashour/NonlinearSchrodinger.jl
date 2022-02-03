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

"""
    function Box(xᵣ::Pair, T; dx, Nₓ, Nₜ, n_periods)
    
Create a `::Box` object with `Nₜ` transverse nodes, `Nₓ` longitudinal nodes (or `dx` grid spacing in the lontidunal direction, only one option can be specified), of size `[(xᵣ.first), (xᵣ.second)]` in the longitudinal direction and `[(-T*n_periods/2), (T*n_periods/2))`
"""
function Box(xᵣ::Pair, T; dx = 0.0, Nₓ = 0, Nₜ = 256, n_periods = 1)
    @info "Initializing Box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ."
    T = n_periods * T
    @info "Longitudinal range is [$(xᵣ.first), $(xᵣ.second)], transverse range is [$(-T/2), $(T/2))"
    dt = T / Nₜ
    t = dt * collect((-Nₜ/2:Nₜ/2-1))
    if dx == 0.0 && Nₓ == 0
        throw(ArgumentError("You must specify either dx or Nₓ"))
    elseif dx != 0.0 && Nₓ != 0
        throw(ArgumentError("You must specify either dx or Nₓ, not both."))
    elseif dx == 0.0 && Nₓ != 0
        dx = (xᵣ.second - xᵣ.first)/Nₓ
        x = collect(xᵣ.first:dx:xᵣ.second)
        Nₓ = length(x)
    elseif dx != 0.0 && Nₓ == 0
        x = collect(xᵣ.first:dx:xᵣ.second)
        Nₓ = length(x)
    end

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
    β::TT
    ψ::Array{Complex{TT}, 2}
    ψ̃::Array{Complex{TT}, 2}
    E::Array{TT, 1}
    PE::Array{TT, 1}
    KE::Array{TT, 1}
    N::Array{TT, 1}
    P::Array{TT, 1}
    save_at::Array{TT, 1}
end # Simulation

"""
    function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}, T̂; α = 0.0, ϵ=0.0, β=0.0) where TT <: Rea
    
Create a `::Sim` object corresponding to eigenvalue `λ` (only used for breathers, ignored for arbitrary initial conditions), initial condition `ψ₀` and algorithm `T̂`. `α` is the Hirota equation parameter and `ϵ` is the Sasa-Satsuma equation parameter. `β` controls pruning.

See also: [`Box`](@ref)
"""
function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}, T̂; α = 0.0, ϵ=0.0, β=0.0, save_interval=:auto) where TT <: Real
    if save_interval==:auto
        save_at = box.x
    else
        si = save_interval*box.dx
        save_at = collect(range(box.x[1],step=si,stop=box.x[end]))
    end
    ψ = Array{Complex{TT}}(undef, box.Nₜ, length(save_at))
    ψ̃ = similar(ψ)
    E = zeros(length(save_at))
    PE = similar(E)    
    KE = similar(E)
    N = similar(E)
    P = similar(E)
    if β < 0.0 
        throw(ArgumentError("αₚ < 0. Set αₚ = 0 to disable pruning or αₚ = Inf for fixed pruning. See documentation)"))
    end
    if β > 0.0 && box.n_periods == 1 
        throw(ArgumentError("Pruning is only applicable when n_periods > 1."))
    end
    # Compute some parameters
    λ, T, Ω = params(λ = λ)
    sim = Sim(λ, T, Ω, box, ψ₀, T̂, α, ϵ, β, ψ, ψ̃, E, PE, KE, N, P, save_at)
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
    f::Dict{Symbol,TT}
    ψ::Matrix{Complex{TT}}
    ψ̃::Matrix{Complex{TT}}
    E::Vector{TT}
    PE::Vector{TT}
    KE::Vector{TT}
    N::Vector{TT}
    P::Vector{TT}
end # Simulation

"""
    function Calc(λ::Array{Complex{TT}}, tₛ, xₛ, seed, box; m=0.0) where TT <: Real 

Create a `::Calc` object with eigenvalues `λ`, shifts `xₛ` and `tₛ` and seeding solution `seed`. `box::Box` is the Box object used for the calculation and `m` is the elliptic parameter for cnoidal solutions. `seed` can have the following values:

`seed = "0"` ``\\implies \\psi_0 = 0``

`seed = "exp"` ``\\implies \\psi_0 = e^{ix}``

`seed = "dn"` ``\\implies \\psi_0 = dn(t, m)e^{ix(1-m/2)}``

`seed = "cn"` ``\\implies \\psi_0 = \\sqrt{m}cn(t, m)e^{ix(m - 1/2)}``

`f::Dict{Symbol, Float64}` is a dictionary of extended NLSE parameters.

See also: [`Box`](@ref)
"""
function Calc(λ::Array{Complex{TT}}, tₛ, xₛ, seed, box; m=0.0, f = Dict(:α=>0.0, :γ=>0.0,:δ=> 0.0)) where TT <: Real
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

    ψ = Array{Complex{TT}}(undef, box.Nₜ, box.Nₓ)
    ψ̃ = similar(ψ)
    E = zeros(box.Nₓ)
    PE = similar(E)    
    KE = similar(E)
    N = similar(E)
    P = similar(E)

    if !haskey(f, :α)
        push!(f, :α => 0.0)
    end
    if !haskey(f, :γ)
        push!(f, :γ => 0.0)
    end
    if !haskey(f, :δ)
        push!(f, :δ => 0.0)
    end

    calc = Calc(λ, T, Ω, χ, tₛ, xₛ, seed, box, m, f, ψ, ψ̃, E, PE, KE, N, P)
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

@memoize function K_cubic(dx::Real, ω)
    ifftshift(cis.(-dx*ω.^2/2))
end
@memoize function K_hirota(dx::Real, ω, α)
    ifftshift(cis.(-dx*(ω.^2/2 .+ α*ω.^3)))
end
@memoize function K_SS(dx::Real, ω, ϵ)
    ifftshift(cis.(-dx*(ω.^2/2 .- ϵ*ω.^3)))
end

function (K̂::Ks)(dx)
    if K̂.α == 0 && K̂.ϵ == 0
        K_cubic(dx, K̂.ω)
    elseif K̂.α != 0 && K̂.ϵ == 0
        K_hirota(dx, K̂.ω, K̂.α)
    elseif K̂.α == 0 && K̂.ϵ != 0
        K_SS(dx, K̂.ω, K̂.ϵ)
    else
        @error "Only ϵ or α can be nonzero at a time"
    end
end

"""
    function Operators(sim)

Construct appropriate operators for `sim::Sim`. These include the Fourier operator `F̂`, inverse Fourier operator `F̃̂`, dispersion operaotr `K̂`, and Hirota/Sasa-Satsuma operators `B̂` and `B̂B`. This object also includes some arrays for temporary operations in Nystrom integrators.

See also: [`Sim`](@ref)
"""
function Operators(sim)
    # Generate FFT Plans to optimize performance
    @info "Generating FFT plans"
    F̂ = plan_fft!(@view sim.ψ[:, 1]) # 26 allocs
    F̃̂ = plan_ifft!(@view sim.ψ[:, 1]) # 34 allocs

    # This stuff should be user controllable
    t_algo = Tsit5()
    t_order = 2
    # Create the integrator for the Burger term
    if sim.α != 0 && sim.ϵ == 0
        D = CenteredDifference(1, t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        function MB!(du, u,p,t)
            du .= 6*sim.α*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(MB!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false, adaptive=false)  
        B̂B = B̂
    elseif sim.α == 0 && sim.ϵ != 0
        D = CenteredDifference(1, t_order, sim.box.dt, sim.box.Nₜ) 
        Q = PeriodicBC(Float64)
        function MS1!(du, u,p,t)
            du .= -6*sim.ϵ*D*Q*u.*abs2.(u)
        end
        prob = ODEProblem(MS1!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false, adaptive=false)  
        function MS2!(du, u,p,t)
            du .= -3*sim.ϵ*D*Q*abs2.(u).*u
        end
        prob = ODEProblem(MS2!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂B = init(prob, t_algo; dt=sim.box.dx,save_everystep=false, adaptive=false)  
    else
        function M0!(du, u,p,t)
            du .= 0*u 
        end
        prob = ODEProblem(M0!, sim.ψ[:, 1], (0, sim.box.dx))
        B̂ = init(prob, t_algo; dt=sim.box.dx,save_everystep=false, adaptive=false)  
        B̂B = B̂
    end
    
    ψ₁ = similar(sim.ψ₀)
    ψ₂ = similar(sim.ψ₀)
    ψ₃ = similar(sim.ψ₀)
    ops = Operators(Ks(sim.α, sim.ϵ, sim.box.ω), F̂, F̃̂, B̂, B̂B, ψ₁, ψ₂, ψ₃)

    return ops
end