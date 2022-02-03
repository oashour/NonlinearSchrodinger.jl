@recipe function f(sim::T, mode=:ψ) where T <: Union{Calc, Sim}
    # Set Defaults
    res_x --> 500
    res_t --> 500
    n_lines --> 10
    skip --> 1

    # Extract Values
    res_x = plotattributes[:res_x]
    res_t = plotattributes[:res_t]
    skip = plotattributes[:skip]
    n_lines = plotattributes[:n_lines]

    # Compute sampling interval
    xₛ = Int(ceil(length(sim.save_at)/res_x))
    tₛ = Int(ceil(sim.box.Nₜ/res_t))


    if mode==:ψ
        # Downsample
        x_ax = sim.box.t[1:tₛ:end]
        y_ax =  sim.save_at[1:xₛ:end] 
        ψ = sim.ψ[1:tₛ:end, 1:xₛ:end]

        xguide --> L"t"
        yguide --> L"x"
        colorbar_title --> L"|\psi|"
        seriescolor --> :viridis
        fontfamily --> "Computer Modern"
        tick_direction --> :out
        s = get(plotattributes, :seriestype, :auto)
        if s == :surface
            zguide --> L"|\psi|"
            colorbar --> false
            camera --> (35, 75)
            grid --> false
        end
        @series begin
            # return series data
            x_ax, y_ax, abs.(ψ)'
        end
    elseif mode==:ψ̃
        s = get(plotattributes, :seriestype, :auto)
        if s == :auto
            # Downsample
            y_ax =  sim.save_at[1:xₛ:end] 
            #ψ̃ = sim.ψ̃[1:tₛ:end, 1:xₛ:end]
            # Extract Lines
            start = sim.box.Nₜ÷2 + 1
            ψ̃ = sim.ψ̃[start:skip:start+skip*n_lines-1, 1:xₛ:end]
            labs = reshape([L"\tilde{\psi}_{%$(i-1)}" for i in 1:skip:skip*n_lines], 1, :)
            tick_direction --> :out
            linewidth --> 1.5
            xguide --> L"x"
            yguide --> L"\log|\tilde{\psi}|"
            fontfamily --> "Computer Modern"
            tick_direction --> :out
            #margin --> 4mm
            legend --> :outertopright
            for i=1:size(ψ̃)[1]
                @series begin
                    label --> labs[i]
                    y_ax,  log.(abs.(ψ̃[i, :]))
                end
            end
        else
            # Downsample
            x_ax = sim.box.ω[1:tₛ:end]
            y_ax =  sim.save_at[1:xₛ:end] 
            ψ̃ = sim.ψ̃[1:tₛ:end, 1:xₛ:end]

            xguide --> L"\omega"
            yguide --> L"x"
            colorbar_title --> L"\log|\tilde{\psi}|"
            seriescolor --> :viridis
            fontfamily --> "Computer Modern"
            tick_direction --> :out
            @series begin
                seriestype := s
                # return series data
                x_ax, y_ax, log.(abs.(ψ̃'))
            end
        end
    elseif mode==:IoM
        # Downscale
        x = sim.save_at[1:xₛ:end]
        E = sim.E[1:xₛ:end]
        KE = sim.KE[1:xₛ:end]
        PE = sim.PE[1:xₛ:end]
        δE = sim.E[1:xₛ:end] .- sim.E[1]
        δP = sim.P[1:xₛ:end] .- sim.P[1]
        δN = sim.N[1:xₛ:end] .- sim.N[1]

        layout --> (2, 2)
        fontfamily --> "Computer Modern"
        tick_direction --> :out
        @series begin
            linewidth --> 1.5
            subplot := 1
            label --> L"E"
            x, E
        end
        @series begin
            linewidth --> 1.5
            subplot := 1
            label --> L"K"
            x, KE
        end
        @series begin
            linewidth --> 1.5
            subplot := 1
            xguide := L"x"
            yguide := L"E"
            label --> L"V"
            x, PE
        end
        @series begin
            linewidth --> 1.5
            subplot := 2
            xguide := L"x"
            yguide := L"\delta E"
            label --> ""
            x, δE
        end
        @series begin
            linewidth --> 1.5
            subplot := 3
            xguide := L"x"
            yguide := L"\delta P"
            label --> ""
            x, δP
        end
        @series begin
            linewidth --> 1.5
            subplot := 4
            xguide := L"x"
            yguide := L"\delta N"
            label --> ""
            x, δN
        end
    else
        @error "Unknown Plotting Mode for NonlinearSchrodinger.jl Objects $mode"
    end

    return nothing
end