@recipe function f(sim::T; x_res=500, t_res=500) where T <: Union{Calc, Sim}
    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    tₛ = Int(ceil(sim.box.Nₜ/t_res))

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

    return nothing
end