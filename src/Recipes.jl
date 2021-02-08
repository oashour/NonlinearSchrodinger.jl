@recipe function f(::Type{T}, sim::T; x_res=500, t_res=500) where T <: Union{Calc, Sim}
    # Compute sampling interval
    xₛ = Int(ceil(sim.box.Nₓ/x_res))
    tₛ = Int(ceil(sim.box.Nₜ/t_res))

    xguide --> L"t"
    yguide --> L"x"
    colorbar_title --> L"|\psi|"
    seriescolor --> :viridis
    #colorbar --> false

    # Downsample
    x_ax = sim.box.t[1:tₛ:end]
    y_ax =  sim.box.x[1:xₛ:end] 
    ψ = sim.ψ[1:tₛ:end, 1:xₛ:end]

    y_yt =  y_ax[1:xₛ:end] 
    yt_loc = Base.OneTo(length(y_ax))


    #yticks --> []
    #yticks --> []
    # Set Axes
    x := x_ax
    y := y_ax
    return abs.(ψ'), x_ax, y_ax
    #return (x_ax, y_ax, abs.(ψ'))
end

#@shorthands density_ψ
#@recipe function f(::Type{Val{:density_ψ}}, plt::AbstractPlot, sim::T; x_res=500, t_res=500) where T
#    # Compute sampling interval
#    xₛ = Int(ceil(sim.box.Nₓ/x_res))
#    tₛ = Int(ceil(sim.box.Nₜ/t_res))
#
#    # Downsample
#    x_ax = sim.box.t[1:tₛ:end]
#    y_ax =  sim.box.x[1:xₛ:end] 
#    ψ = sim.ψ[1:tₛ:end, 1:xₛ:end]
#    @series begin
#        seriestype := :heatmap
#        x, y, z
#    end
#end