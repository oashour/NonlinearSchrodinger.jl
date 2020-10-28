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
    if sim.algorithm.variant === :A
        soln_loop_A(sim, ops, ind_p)
    elseif sim.algorithm.variant === :B
        soln_loop_B(sim, ops, ind_p)
    end

    @info "Computation Done!"
end #solve

function soln_loop_A(sim, ops, ind_p)
    F̂ = plan_fft(@view sim.ψ[:, 1]) # 26 allocs
    @progress for i = 1:sim.box.Nₓ-1
        @time @views sim.ψ[:, i+1] .= sim.algorithm.T̂(sim.ψ[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            ops.F̂*view(sim.ψ,:,i+1)
            for j in ind_p
                sim.ψ[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ[j, i+1]))
            end
            ops.F̃̂*view(sim.ψ,:,i+1)
        end # if
        @views sim.ψ̃[:, i+1] .= fftshift(F̂*sim.ψ[:, i+1])/sim.box.Nₜ
    end # for
end

function soln_loop_B(sim, ops, ind_p)
    F̃̂ = plan_ifft(@view sim.ψ[:, 1]) # 34 allocs
    F̂ = plan_fft(@view sim.ψ[:, 1]) # 26 allocs
    sim.ψ̃[:, 1] .= F̂*view(sim.ψ, :, 1)
    @progress for i = 1:sim.box.Nₓ-1
        @time sim.ψ̃[:, i+1] .= sim.algorithm.T̂(sim.ψ̃[:, i], sim.box.dx, ops)
        # Pruning
        if sim.αₚ > 0 
            for j in ind_p
                sim.ψ̃[j, i+1] *= exp(-sim.αₚ*abs(sim.ψ̃[j, i+1]))
            end
        end 
        sim.ψ[:, i+1] .= F̃̂*view(sim.ψ̃, :, i+1)
    end
    sim.ψ̃ .= fftshift(sim.ψ̃, 1)
end