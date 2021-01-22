include("CubicSolvers.jl")
include("HirotaSolvers.jl")
include("SSSolvers.jl")

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
    ind_p = zeros(Int64, 1)
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
        @debug "Computed pruning indices $ind_p"
    end
    
    ops = Operators(sim)

    @info "Starting evolution"
    # Define these tuples as global consts
    if sim.T̂ ∈ (T1A!, T2A!, T4A_TJ!, T6A_TJ!, T8A_TJ!, T4A_SF!, T4A_SF!, T6A_SF!, T8A_SF!, T4A_N!, T6A_N!, T8A_N!, T6A_OP!, T8A_OP!, T1A_H!, T2A_H!, T1A_SS!)
        soln_loop_A(sim, ops, ind_p)
    elseif sim.T̂ ∈ (T1B!, T2B!, T4B_TJ!, T6B_TJ!, T8B_TJ!, T4B_SF!, T6B_SF!, T8B_SF!, T4B_N!, T6B_N!, T8B_N!, T6B_OP!, T8B_OP!)
        soln_loop_B(sim, ops, ind_p)
    else
        throw(ArgumentError("Unknown algorithm $(sim.T̂)"))
    end

    @info "Computation Done!"
end #solve

function soln_loop_A(sim, ops, ind_p)
    @progress for i = 1:sim.box.Nₓ-1
        @views sim.T̂(sim.ψ[:, i+1], sim.ψ[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            ops.F̂*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ[j, i+1]))
            end
            ops.F̃̂*view(sim.ψ,:,i+1)
        end # if
    end # for
    @info "Computing Spectrum"
    sim.ψ̃ .= fftshift(fft(sim.ψ, 1), 1)./sim.box.Nₜ
end

function soln_loop_B(sim, ops, ind_p)
    F̃̂ = plan_ifft(@view sim.ψ[:, 1]) 
    F̂ = plan_fft(@view sim.ψ[:, 1])
    sim.ψ̃[:, 1] .= F̂*view(sim.ψ, :, 1)
    @progress for i = 1:sim.box.Nₓ-1
        @views sim.T̂(sim.ψ̃[:, i+1], sim.ψ̃[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            for j in ind_p
                sim.ψ̃[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ̃[j, i+1]))
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