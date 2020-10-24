module Utilities
using FFTW, JLD
export compute_energy!, compute_spectrum!

"""
    function compute_spectrum!(obj)

Computes the normalized spectrum of `obj.ψ` with the center frequency shifted to the center
and saves it in `obj.ψ̃`. `obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_ψ̃`](@ref)
"""
function compute_spectrum!(obj)
    println("==========================================")
    println("Computing spectrum")
    if obj.solved
        obj.ψ̃ = fftshift(fft(obj.ψ, 1), 1)/obj.box.Nₜ
        obj.spectrum_computed = true
    else
        throw(ArgumentError("Trying to compute spectrum of an unsolved simulation. Please solve the model first."))
    end
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
    if ~obj.spectrum_computed
        println("IoM calculation requested without a spectrum calculation.")
        compute_spectrum!(obj)
    end
    println("==========================================")
    println("Computing integrals of motions")
    # Compute the energies
    # Do I need to find a better way of doing these integrals?
    # We should have 
    # sim.norm = sum(abs.(sim.ψ).^2, dims=2)[:]*sim.box.dt/sim.box.T but dt/T = 1/Nt, thus
    obj.N = sum(abs2.(obj.ψ), dims=1)[:]/obj.box.Nₜ
    obj.PE = -0.5*sum(abs2.(obj.ψ).^2,dims=1)[:]./(obj.N*obj.box.Nₜ)
    obj.KE = 0.5*sum((obj.box.ω.^2) .* (abs2.(obj.ψ̃)),dims=1)[:]./obj.N
    obj.P = -imag.(sum(im * (obj.box.ω) .* (abs2.(obj.ψ̃)),dims=1)[:]./obj.N)
    obj.E = obj.KE + obj.PE
    obj.dE = obj.E .- obj.E[1]
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