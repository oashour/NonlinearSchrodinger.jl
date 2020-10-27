include("CubicSolvers.jl")
include("HirotaSolvers.jl")

struct Sim{TT<:Real}
    λ::Complex{TT}
    T::TT
    Ω::TT
    box::Box{TT}
    ψ₀::Array{Complex{TT}, 1}
    t_order::Int64
    x_order::Int64
    α::TT
    αₚ::TT
    ψ::Array{Complex{TT}, 2}
    ψ̃::Array{Complex{TT}, 2}
    E::Array{TT, 1}
    PE::Array{TT, 1}
    KE::Array{TT, 1}
    #dE::Array{TT, 1}
    N::Array{TT, 1}
    P::Array{TT, 1}
end # Simulation

function Sim(λ, box::Box, ψ₀::Array{Complex{TT}, 1}; t_order = 2, x_order = 2, α = 0.0, αₚ = 0.0) where TT <: Real
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
    sim = Sim(λ, T, Ω, box, ψ₀, t_order, x_order, α, αₚ, ψ, ψ̃, E, PE, KE, N, P)
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
    if sim.x_order == 1 && sim.α == 0
        T̂ = T₁ˢ 
    elseif sim.x_order == 2 && sim.α == 0
        T̂ = T₂ˢ
    elseif sim.x_order == 4 && sim.α == 0
        T̂ = T₄ˢ 
    elseif sim.x_order == 6 && sim.α == 0
        T̂ = T₆ˢ 
    elseif sim.x_order == 8 && sim.α == 0
        T̂ = T₈ˢ 
    elseif sim.x_order == 1 && sim.α >= 0
        T̂ = T₁ʰ 
    elseif sim.x_order == 2 && sim.α >= 0
        T̂ = T₂ʰ 
    elseif sim.x_order == 4 && sim.α >= 0
        T̂ = T₄ʰ 
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

"""
    solve!(sim::Simulation)

Solves the `Simulation` object `sim` using the techniques its attributes specify.

See also: [`init_sim`](@ref), [`NLSS.Plotter.plot_ψ`](@ref)
"""
function solve!(sim::Sim)
    @info println("Solving cubic NLSE with the following options:")
    print(sim)       

    # Find
    sim.ψ[:, 1] = sim.ψ₀
    # Check for pruning and calculate indices
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
        @debug "Computed pruning indices $ind_p"
    end
    
    ops = Operators(sim)

    @info "Starting evolution"
    soln_loop(sim, ops)
    #sim.solved = true

    @info "Computation Done!"

    return nothing
end #solve

function soln_loop(sim, ops)
    @progress for i = 1:sim.box.Nₓ-1
        @views sim.ψ[:, i+1] .= ops.T̂(sim.ψ[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            ops.F̂*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ[j, i+1]))
            end
            F̃̂*view(sim.ψ,:,i+1)
        end # if
    end # for
end