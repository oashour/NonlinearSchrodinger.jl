module Plotter
using ..Simulation
using Plots; gr()
Plots.GRBackend()
using LaTeXStrings
using FFTW

export plot_ψ, plot_ψ̃, plot_CoM

"""
    function plot_CoM(sim::Sim, x_res::Int64 = 500)

Plots the integrals of motion for a `Sim` object `sim` with a resolution `x_res` points in
the x-direction. Produces plots of the energy, kinetic energy, potential energy, energy
error, particle number and momentum

See also: [`Simulation.compute_CoM!`](@ref)
"""
function plot_CoM(sim::Sim, x_res = 500)
    println("Plotting energy with a resolution of $x_res")
    xₛ = Int(ceil(sim.box.Nₜ/x_res))
    p1 = plot(sim.box.x[1:xₛ:end], [sim.E[1:xₛ:end], sim.KE[1:xₛ:end], sim.PE[1:xₛ:end]], label = [L"E" L"T" L"V"])
    xlabel!(L"x")
    ylabel!(L"E")
    xlims!((minimum(sim.box.x), maximum(sim.box.x)))
    p2 = plot(sim.box.x[1:xₛ:end], sim.dE[1:xₛ:end], label = "")
    xlabel!(L"x")
    ylabel!(L"\delta E")
    xlims!((minimum(sim.box.x), maximum(sim.box.x)))
    p3 = plot(sim.box.x[1:xₛ:end], [sim.P[1:xₛ:end]], label = "")
    xlabel!(L"x")
    ylabel!(L"P")
    xlims!((minimum(sim.box.x), maximum(sim.box.x)))
    p4 = plot(sim.box.x[1:xₛ:end], [sim.N[1:xₛ:end]], label = "")
    xlabel!(L"x")
    ylabel!(L"N")
    xlims!((minimum(sim.box.x), maximum(sim.box.x)))

    return p1, p2, p3, p4
end

"""
    plot_ψ(sim::Sim; mode = "density", power=1, x_res=500, t_res=512)

Plot `sim.abs.(ψ).^power` on a grid of size approximately `t_res`x`x_res`
The function can plot either a heatmap (`mode = "density"`) or a 3D surface
(`mode = "surface"`)

**TODO:** Insert pictures here

See also: [`Simulation.solve!`](@ref)
"""
function plot_ψ(sim::Sim; mode = "density", power=1, x_res=500, t_res=512)
    if ~sim.solved
        throw(ArgumentError("The simulation has not been solved, unable to plot."))
    end
    println("==========================================")
    println("Plotting |ψ|^$power in $mode mode with ~$x_res longitudinal and ~$t_res transverse points.")

    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    tₛ = Int(ceil(sim.box.Nₜ/t_res))
    # Downsample
    t = sim.box.t[1:tₛ:end]
    x = sim.box.x[1:xₛ:end]
    ψ = sim.ψ[1:xₛ:end, 1:tₛ:end]

    # Figure out labels
    if power == 1
        clabel = L"|\psi|"
    elseif power == 2
        clabel = L"|\psi|^2"
    end

    # Plot
    if mode == "density"
        p = heatmap(t, x, abs.(ψ).^power, tick_direction=:out, colorbar_title=clabel)
    elseif mode == "surface"
        p = surface(t, x, abs.(ψ).^power, zlabel=clabel, camera=(35, 68))
    end
    # Adjust Attributes
    xlabel!(L"x")
    ylabel!(L"t")
    xlims!((-sim.params.T/2*sim.box.n_periods, sim.params.T/2*sim.box.n_periods))
    ylims!((minimum(sim.box.x), maximum(sim.box.x)))

    println("Plotting done!")
    println("==========================================")

    return p
end #plot_ψ

"""
    plot_ψ̃(sim::Sim; mode = "density", power=1, x_res=500, t_res=512)

Plot `log(sim.abs.(ψ))` on a grid of size approximately `ω_res`x`x_res`
if `mode = density`, or with `n_lines` total, skipping `skip` between each 
line if `mode = lines`. The latter uses `x_res` points as well.

**TODO:** Insert pictures here

See also: [`Simulation.compute_spectrum!`](@ref)
"""
function plot_ψ̃(sim::Sim; mode = "density", x_res=500, ω_res=512, skip = 1, n_lines = 10)
    if ~sim.solved
        throw(ArgumentError("The simulation has not been solved, unable to plot."))
    end
    println("==========================================")
    println("Plotting log(|ψ̃|) in $mode mode with ~$x_res longitudinal and ~$ω_res transverse points.")

    # Figure out labels
    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    ωₛ = Int(ceil(sim.box.Nₜ/ω_res))
    # Plot
    if mode == "density"
        # Compute sampling interval
        xₛ = Int(ceil(sim.box.Nₓ/x_res))
        ωₛ = Int(ceil(sim.box.Nₜ/ω_res))
        # Downsample
        ω = sim.box.ω[1:ωₛ:end]
        x = sim.box.x[1:xₛ:end]
        ψ̃ = sim.ψ̃[1:xₛ:end, 1:ωₛ:end]

        p = heatmap(ω, x, log.(abs.(ψ̃)), tick_direction=:out, colorbar_title=L"\log|\tilde{\psi}|")
        xlabel!(L"\omega")
        ylabel!(L"x")
        xlims!((-sim.box.Nₜ/2*sim.params.Ω/sim.box.n_periods, sim.box.Nₜ/2*sim.params.Ω/sim.box.n_periods))
        ylims!((minimum(sim.box.x), maximum(sim.box.x)))
    elseif mode == "lines"
        ψ̃ = ifftshift(sim.ψ̃, 2)[1:xₛ:end, 1:skip:skip*n_lines]
        x = sim.box.x[1:xₛ:end]
        labs = reshape([L"\tilde{\psi}_{%$(i-1)}" for i in 1:skip:skip*n_lines], 1, :)
        p = plot(x, log.(abs.(ψ̃)), label = labs)
        xlabel!(L"x")
        ylabel!(L"\log|\tilde{\psi}|")
    end

    println("Plotting done!")
    println("==========================================")

    return p
end #plot_ψ̃
end