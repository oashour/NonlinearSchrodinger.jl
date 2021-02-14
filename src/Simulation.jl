include("CubicSolvers.jl")
include("HirotaSolvers.jl")
include("SSSolvers.jl")

"""
    solve!(sim::Sim)

Solves the object `sim::Sim`.

See also: [`Sim`](@ref)
"""
function solve!(sim::Sim)
    # Set initial condition
    sim.ψ[:, 1] = sim.ψ₀
    # Check for pruning and calculate indices
    ind_p = zeros(Int64, 1)
    if sim.β > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
        @debug "Computed pruning indices $ind_p"
    end
    
    ops = Operators(sim)

    @info "Starting evolution"
    # Define these tuples as global consts
    if sim.T̂ ∈ (T1A!, T2A!, T4A_TJ!, T6A_TJ!, T8A_TJ!, T4A_SF!, T4A_SF!, T6A_SF!, T8A_SF!, T4A_CMP!, T6A_CMP!, T8A_CMP!, T6A_Ss14!, T6A_Ys7!, T6A_KLs9!, T8A_Ss15!, T8A_KLs17!, T1A_H!, T2A_H!, T1A_SS!)
        soln_loop_A(sim, ops, ind_p)
    elseif sim.T̂ ∈ (T1B!, T2B!, T4B_TJ!, T6B_TJ!, T8B_TJ!, T4B_SF!, T6B_SF!, T8B_SF!, T4B_CMP!, T6B_CMP!, T8B_CMP!, T6B_Ss14!, T6B_Ys7!, T6B_KLs9!, T8B_Ss15!, T8B_KLs17!)
        soln_loop_B(sim, ops, ind_p)
    else
        throw(ArgumentError("Unknown algorithm $(sim.T̂)"))
    end

    @info "Computation Done!"
end #solve

"""
    soln_loop_A(sim, ops, ind_p))

Propagate in `x` for integrators of type `A` (i.e. starting with nonlinear step, thus in real space).

See also: [`Sim`](@ref), [`solve!`](@ref)
"""
function soln_loop_A(sim, ops, ind_p)
    for i = 1:sim.box.Nₓ-1
        @views sim.T̂(sim.ψ[:, i+1], sim.ψ[:, i], sim.box.dx, ops)
        # Pruning
        if sim.β > 0 
            ops.F̂*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.β*abs(sim.ψ[j, i+1]))
            end
            ops.F̃̂*view(sim.ψ,:,i+1)
        end # if
    end # for
    @info "Computing Spectrum"
    sim.ψ̃ .= fftshift(fft(sim.ψ, 1), 1)./sim.box.Nₜ
end

"""
    soln_loop_B(sim, ops, ind_p))

Propagate in `x` for integrators of type `B` (i.e. starting with dispersion step, thus in frequency space).

See also: [`Sim`](@ref), [`solve!`](@ref)
"""
function soln_loop_B(sim, ops, ind_p)
    F̃̂ = plan_ifft(@view sim.ψ[:, 1]) 
    F̂ = plan_fft(@view sim.ψ[:, 1])
    sim.ψ̃[:, 1] .= F̂*view(sim.ψ, :, 1)
    for i = 1:sim.box.Nₓ-1
        @views sim.T̂(sim.ψ̃[:, i+1], sim.ψ̃[:, i], sim.box.dx, ops)
        # Pruning
        if sim.β > 0 
            for j in ind_p
                sim.ψ̃[j, i+1] *= exp(-sim.β*abs(sim.ψ̃[j, i+1]))
            end
        end 
        # Compute the actual ψ
        #sim.ψ[:, i+1] .= F̃̂*view(sim.ψ̃, :, i+1)
    end
    # Shift the FT
    @info "Computing Spectrum"
    sim.ψ .= ifft(sim.ψ̃, 1)
    sim.ψ̃ .= fftshift(sim.ψ̃, 1)./sim.box.Nₜ
end