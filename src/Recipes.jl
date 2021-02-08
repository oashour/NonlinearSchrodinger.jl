@recipe function f(sim::T, mode=:ψ; x_res=500, t_res=500, n_lines=10, skip=1) where T <: Union{Calc, Sim}
    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    tₛ = Int(ceil(sim.box.Nₜ/t_res))

    if mode==:ψ
        # Downsample
        x_ax = sim.box.t[1:tₛ:end]
        y_ax =  sim.box.x[1:xₛ:end] 
        ψ = sim.ψ[1:tₛ:end, 1:xₛ:end]

        xguide --> L"t"
        yguide --> L"x"
        colorbar_title --> L"|\psi|"
        seriescolor --> :viridis
        @series begin
            # return series data
            x_ax, y_ax, abs.(ψ)'
        end
    elseif mode==:ψ̃
        # Downsample
        x_ax = sim.box.ω[1:tₛ:end]
        y_ax =  sim.box.x[1:xₛ:end] 
        ψ̃ = sim.ψ̃[1:tₛ:end, 1:xₛ:end]

        xguide --> L"\omega"
        yguide --> L"x"
        colorbar_title --> L"\log|\tilde{\psi}|"
        seriescolor --> :viridis
        @series begin
            # return series data
            x_ax, y_ax, log.(abs.(ψ̃'))
        end
    elseif mode==:ψ̃_lines
        # Downsample
        y_ax =  sim.box.x[1:xₛ:end] 
        #ψ̃ = sim.ψ̃[1:tₛ:end, 1:xₛ:end]
        # Extract Lines
        start = sim.box.Nₜ÷2 + 1
        ψ̃ = sim.ψ̃[start:skip:start+skip*n_lines-1, 1:xₛ:end]
        labs = reshape([L"\tilde{\psi}_{%$(i-1)}" for i in 1:skip:skip*n_lines], 1, :)
        tick_direction --> :out
        linewidth --> 1.5
        xguide --> L"x"
        yguide --> L"\log|\tilde{\psi}|"
        margin --> 4mm
        legend --> :outertopright
        for i=1:size(ψ̃)[1]
            @series begin
                label --> labs[i]
                y_ax,  log.(abs.(ψ̃[i, :]))
            end
        end
    end
    return nothing
end