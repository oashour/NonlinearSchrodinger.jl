include("CubicSolvers.jl")
include("HirotaSolvers.jl")

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
    ind_p = []
    if sim.αₚ > 0 
        ind_p = [i for i in 2:(sim.box.Nₜ÷2+1) if (i-1)%sim.box.n_periods != 0]
        ind_p = sort([ind_p; sim.box.Nₜ.-ind_p.+2])
        @debug "Computed pruning indices $ind_p"
    end
    
    ops = Operators(sim)

    @info "Starting evolution"
    if sim.variant == "A" || sim.variant == "FA"
        soln_loop_A(sim, ops, ind_p)
    elseif sim.variant == "B"
        soln_loop_B(sim, ops, ind_p)
    end
   #sim.solved = true

    @info "Computation Done!"

    return nothing
end #solve

function soln_loop_A(sim, ops, ind_p)
    @progress for i = 1:sim.box.Nₓ-1
        @views sim.ψ[:, i+1] .= ops.T̂(sim.ψ[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            ops.F̂*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ[j, i+1]))
            end
            ops.F̃̂*view(sim.ψ,:,i+1)
        end # if
    end # for
end

function soln_loop_B(sim, ops, ind_p)
    F̂ = plan_fft(@view sim.ψ[:, 1]) # 26 allocs
    F̃̂ = plan_ifft(@view sim.ψ[:, 1]) # 34 allocs
    sim.ψ̃[:, 1] .= F̂*view(sim.ψ, :, 1)
    @progress for i = 1:sim.box.Nₓ-1
        @views sim.ψ̃[:, i+1] .= ops.T̂(sim.ψ̃[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            for j in ind_p
                sim.ψ̃[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ̃[j, i+1]))
            end
        end # if
        sim.ψ[:, i+1] .= F̃̂*view(sim.ψ̃, :, i+1)
    end # for
end