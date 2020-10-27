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
    # Needs optimization
    obj.ψ̃ .= fftshift(fft(obj.ψ, 1), 1)/obj.box.Nₜ
    println("Spectrum computed")
    println("==========================================")
end

"""
    function compute_IoM!(obj)

Computes the integrals of motion of `obj.ψ` and saves them in respective fields of `obj`.
`obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_CoM`](@ref)
"""
function compute_IoM!(obj)
    println("==========================================")
    println("Computing integrals of motions")
    # Compute the energies
    # Do I need to find a better way of doing these integrals?
    # We should have 
    # sim.norm = sum(abs.(sim.ψ).^2, dims=2)[:]*sim.box.dt/sim.box.T but dt/T = 1/Nt, thus
    obj.N .= sum(abs2.(obj.ψ), dims=1)[:]/obj.box.Nₜ
    obj.PE .= -0.5*sum(abs2.(obj.ψ).^2,dims=1)[:]./(obj.N*obj.box.Nₜ)
    obj.KE .= 0.5*sum((obj.box.ω.^2) .* (abs2.(obj.ψ̃)),dims=1)[:]./obj.N
    obj.P .= -imag.(sum(im * (obj.box.ω) .* (abs2.(obj.ψ̃)),dims=1)[:]./obj.N)
    obj.E .= obj.KE + obj.PE
    #obj.dE .= obj.E .- obj.E[1]
    println("Integrals of motion computed.")
    println("==========================================")
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