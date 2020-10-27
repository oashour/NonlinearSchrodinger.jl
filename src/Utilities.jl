module Utilities
using FFTW, JLD
export compute_energy!, compute_spectrum!, Box

struct Box{TT<:Real}
    t::Array{TT, 1}
    ω::Array{TT, 1}
    x::Array{TT, 1}
    Nₜ::Int64
    Nₓ::Int64
    dt::TT
    dx::TT
    n_periods::Int64
end #SimulationBox

function Box(xᵣ::Pair, T; dx = 1e-3, Nₜ = 256, n_periods = 1)
    println("==========================================")
    println("Initializing simulation box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ.")
    T = n_periods * T
    println("Longitudinal range is [$(xᵣ.first), $(xᵣ.second)], transverse range is [$(-T/2), $(T/2))")
    dt = T / Nₜ
    t = dt * collect((-Nₜ/2:Nₜ/2-1))

    x = collect(xᵣ.first:dx:xᵣ.second)
    Nₓ = length(x)
    ω = 2π/T * collect((-Nₜ/2:Nₜ/2-1))

    box = Box(t, ω, x, Nₜ, Nₓ, dt, dx, n_periods)

    println("Done computing t, x, ω")
    println("==========================================")

    return box
end

"""
    function compute_spectrum!(obj)

Computes the normalized spectrum of `obj.ψ` with the center frequency shifted to the center
and saves it in `obj.ψ̃`. `obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_ψ̃`](@ref)
"""
function compute_spectrum!(obj)
    println("==========================================")
    println("Computing spectrum")

    FFTW.set_num_threads(4)
    F̂ = plan_fft(obj.ψ, 1) # 26 allocs
    obj.ψ̃ .= fftshift(F̂*obj.ψ, 1)/obj.box.Nₜ
    FFTW.set_num_threads(1) # set it back to 1, not useful in smaller FFTs

    println("Spectrum computed")
    println("==========================================")
end

"""
    function compute_IoM!(obj)

Computes the integrals of motion of `obj.ψ` and saves them in respective fields of `obj`.
`obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_CoM`](@ref)
"""
function compute_IoM!(obj; normalize=false)
    #println("==========================================")
    #println("Computing integrals of motions")
    ψ² = abs2.(obj.ψ)
    ψ̃² = abs2.(obj.ψ̃)
    # The norm is preserved extremely well so if the W.F. is normalized we do not worry about dividing the IoM by the norm to improve performance
    if normalize
        obj.N .= sum(ψ², dims=1)[:]./obj.box.Nₜ
        obj.PE .= -0.5*sum(ψ².^2,dims=1)[:]./(obj.N*obj.box.Nₜ)
        obj.KE .= 0.5*sum((obj.box.ω.^2 .* ψ̃²),dims=1)[:]./obj.N
        obj.P .= imag.(sum(im * obj.box.ω.* ψ̃²,dims=1)[:]./obj.N)
        obj.E .= obj.KE .+ obj.PE
    else
        obj.N .= sum(ψ², dims=1)[:]./obj.box.Nₜ
        obj.PE .= -0.5*sum(ψ².^2,dims=1)[:]./obj.box.Nₜ
        obj.KE .= 0.5*sum((obj.box.ω.^2 .* ψ̃²),dims=1)[:]
        obj.P .= imag.(sum(im * obj.box.ω.* ψ̃²,dims=1)[:])
        obj.E .= obj.KE .+ obj.PE
    end

    #println("Integrals of motion computed.")
    #rintln("==========================================")
end

function save(obj, filename)
    filename = string(filename, ".jld")
    println(filename)
    jldopen(filename, "w") do file
        write(file, "result", obj)
    end
    println("Saved to file $filename")
    
    return nothing
end

function load(filename)
    filename = string(filename, ".jld")
    println(filename)
    obj = jldopen(filename, "r") do file
        read(file, "result")
    end
    println("Loaded from file $filename")

    return obj
end

end #module