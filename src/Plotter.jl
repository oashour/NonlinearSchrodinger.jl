module Plotter
using ..Simulation
using Plots; gr()
using Plots.PlotMeasures
using LaTeXStrings
using FFTW

export plot_ψ, plot_ψ̃, plot_CoM

"""
    function plot_IoM(sim::Sim, x_res::Int64 = 500)

Plots the integrals of motion for a `Sim` object `sim` with a resolution `x_res` points in
the x-direction. Produces plots of the energy, kinetic energy, potential energy, energy
error, particle number and momentum

See also: [`Simulation.compute_CoM!`](@ref)
"""
function plot_IoM(obj; x_res = 500, mode = "merged")
    @info "Plotting IoM with a resolution of $x_res"
    xₛ = Int(ceil(obj.box.Nₓ/x_res))
    x = obj.box.x[1:xₛ:end]
    p1 = plot(x, [obj.E[1:xₛ:end], obj.KE[1:xₛ:end], obj.PE[1:xₛ:end]], label = [L"E" L"T" L"V"])
    xlabel!(L"x")
    ylabel!(L"E")
    xlims!((minimum(obj.box.x), maximum(obj.box.x)))
    δE = obj.E[1:xₛ:end] .- obj.E[1]
    p2 = plot(x, δE, label = "")
    xlabel!(L"x")
    ylabel!(L"\delta E")
    xlims!((minimum(obj.box.x), maximum(obj.box.x)))
    δP = obj.P[1:xₛ:end] .- obj.P[1]
    p3 = plot(x, δP, label = "")
    xlabel!(L"x")
    ylabel!(L"\delta P")
    xlims!((minimum(obj.box.x), maximum(obj.box.x)))
    δN = obj.N[1:xₛ:end] .- obj.N[1]
    p4 = plot(x, δN, label = "")
    xlabel!(L"x")
    ylabel!(L"\delta N")
    xlims!((minimum(obj.box.x), maximum(obj.box.x)))

    if mode == "merged"
        p = plot(p1, p2, p3, p4, layout=(2,2), link = :x)
    elseif mode == "separate"
        p = (p1, p2, p3, p4)
    end

    @info "Plotting IoM Done"
    return p
end

"""
    plot_ψ(sim::Sim; mode = "density", power=1, x_res=500, t_res=512)

Plot `sim.abs.(ψ).^power` on a grid of size approximately `t_res`x`x_res`
The function can plot either a heatmap (`mode = "density"`) or a 3D surface
(`mode = "surface"`)

**TODO:** Insert pictures here

See also: [`Simulation.solve!`](@ref)
"""
function plot_ψ(sim; mode = "density", power=1, x_res=500, t_res=512)
    @info "Plotting |ψ|^$power in $mode mode with ~$x_res longitudinal and ~$t_res transverse points."

    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    tₛ = Int(ceil(sim.box.Nₜ/t_res))
    # Downsample
    t = sim.box.t[1:tₛ:end]
    x = sim.box.x[1:xₛ:end]
    ψ = sim.ψ[1:tₛ:end, 1:xₛ:end]

    # Figure out labels
    if power == 1
        clabel = L"|\psi|"
    elseif power == 2
        clabel = L"|\psi|^2"
    end

    # Plot
    if mode == "density"
        p = heatmap(t, x, abs.(ψ').^power, 
                    tick_direction=:out, 
                    colorbar_title=clabel,
                    margin = 4mm,
                    xlabel=L"t",
                    ylabel=L"x")
    elseif mode == "surface"
        p = surface(t, x, abs.(ψ').^power, 
                    colorbar_title=clabel, 
                    camera=(20, 68),
                    grid=false,
                    xlabel=L"t",
                    ylabel=L"x")
    elseif mode == "contourf"
        # Should also adjust colorbar ticks when the feature is added
        p = contourf(t, x, abs.(ψ').^power, 
                    tick_direction=:out, 
                    colorbar_title=clabel,
                    margin = 4mm,
                    xlabel=L"t",
                    ylabel=L"x",
                    levels=[0.25,0.4,0.8,1,1.25,1.5,1.75],
                    clims=(0.25,1.75))
    elseif mode == "contourf_clean"
        p = contourf(t, x, abs.(ψ').^power, 
                    showaxis=false,
                    ticks=[],
                    levels=[0.2,0.4,0.8,1,1.25,1.5,1.8],
                    colorbar=false)
    elseif mode == "contourf_inset"
        l = @layout [
            a 
            d{0.793h}
        ]
        ps = surface(t, x, abs.(ψ').^power, 
                    colorbar_title=clabel,
                    zlabel=clabel,
                    camera=(20, 68),
                    grid=false,
                    xlabel=L"t",
                    ylabel=L"x")
        pe = contourf(0,
                    showaxis=false,
                    colorbar=false,
                    ticks=[],
                    grid=false)
        p = plot(pe, ps, layout = l, size=(600, 500), margins=1mm)
        pc = contourf!(t, x, abs.(ψ').^power, 
                    inset = (2, bbox(-0.09, -0.39, 0.55, 0.55)), 
                    subplot=3,
                    showaxis=false,
                    ticks=[],
                    levels=[0.2,0.4,0.8,1,1.25,1.5,1.8],
                    colorbar=false)

    elseif mode == "contourf_top"
        l = @layout [
             [a{0.01w} b{0.89w} c{0.1w}]
            d{0.6666667h}
        ]
        ps = surface(t, x, abs.(ψ').^power, 
                    colorbar_title=clabel,
                    zlabel=clabel,
                    camera=(20, 68),
                    grid=false,
                    xlabel=L"t",
                    ylabel=L"x")
        pe1 = contourf(0,
                    showaxis=false,
                    colorbar=false,
                    ticks=[],
                    grid=false)
        pc = contourf(t, x, abs.(ψ').^power, 
                    showaxis=false,
                    ticks=[],
                    levels=[0.2,0.4,0.8,1,1.25,1.5,1.8],
                    colorbar=false)
        pe2 = contourf(0,
                    showaxis=false,
                    colorbar=false,
                    ticks=[],
                    grid=false)

        p = plot(pe1, pc, pe2, ps, layout = l, size=(600, 600))
    end
    # Adjust Attributes
    xlims!((-sim.T/2*sim.box.n_periods, sim.T/2*sim.box.n_periods))
    ylims!((minimum(sim.box.x), maximum(sim.box.x)))

    @info "Plotting ψ done!"
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
function plot_ψ̃(sim; mode = "density", x_res=500, ω_res=512, skip = 1, n_lines = 10)
    @info "Plotting log(|ψ̃|) in $mode mode with ~$x_res longitudinal and ~$ω_res transverse points."

    # Figure out labels
    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    ωₛ = Int(ceil(sim.box.Nₜ/ω_res))
    n = Int.(round.(sim.box.ω/sim.Ω)) 
    n = n[1:ωₛ:end]
    ω = sim.box.ω[1:ωₛ:end]
    x = sim.box.x[1:xₛ:end]
    # Plot
    if mode == "density"
        # Downsample
        ψ̃ = sim.ψ̃[1:ωₛ:end, 1:xₛ:end]

        #heatmap(n, x, log.(abs.(ψ̃')), xmirror = true)
        p = heatmap(ω, x, log.(abs.(ψ̃')), 
                    tick_direction=:out, 
                    colorbar_title=L"\log|\tilde{\psi}|",
                    xlabel = L"\omega",
                    ylabel = L"x",
                    margin = 4mm)
        xlims!((-sim.box.Nₜ/2*sim.Ω/sim.box.n_periods, sim.box.Nₜ/2*sim.Ω/sim.box.n_periods))
        ylims!((minimum(sim.box.x), maximum(sim.box.x)))
    elseif mode == "lines"
        start = sim.box.Nₜ÷2 + 1
        ψ̃ = sim.ψ̃[start:skip:start+skip*n_lines-1, 1:xₛ:end]
        labs = reshape([L"\tilde{\psi}_{%$(i-1)}" for i in 1:skip:skip*n_lines], 1, :)
        p = plot(x, log.(abs.(ψ̃')), 
                    label = labs,  
                    legend = :outertopright,
                    tick_direction = :out,
                    linewidth=1.5,
                    xlabel = L"x",
                    ylabel = L"\log|\tilde{\psi}|",
                    margin = 4mm)
        xlims!((minimum(sim.box.x), maximum(sim.box.x)))
    end

    @info "Plotting ψ̃ done!"

    return p
end #plot_ψ̃
end