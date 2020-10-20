export print, ψ₀_periodic, compute_energy!, compute_spectrum!

"""
    function print(sim::Sim)

Prints information about the `Sim` instance `sim`
"""
function print(sim::Sim)
    println("Box Properties:")
    println("------------------------------------------")
    println("dx = $(sim.box.dx) (Nₓ = $(sim.box.Nₓ))")
    println("Nₜ = $(sim.box.dx) (dt = $(sim.box.dt))")
    println("Parameters:")
    println("------------------------------------------")
    println("λ = $(sim.λ)")
    println("Ω = $(sim.Ω)")
    println("T = $(sim.T)")
    println("------------------------------------------")
    # Should add information about ψ₀
    if sim.algorithm == "2S"
        println("Algorithm: second order symplectic")
    elseif sim.algorithm == "4S"
        println("Algorithm: fourth order symplectic")
    elseif sim.algorithm == "6S"
        println("Algorithm: sixth order symplectic")
    elseif sim.algorithm == "8S"
        println("Algorithm: eighth order symplectic")
    else
        throw(ArgumentError("Algorithm type unknown, please check the documentation"))
    end
    println("------------------------------------------")
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
    if obj.solved
        obj.ψ̃ = fftshift(fft(obj.ψ, 2), 2)/obj.box.Nₜ
        obj.spectrum_computed = true
    else
        throw(ArgumentError("Trying to compute spectrum of an unsolved simulation. Please solve the model first."))
    end
    println("Spectrum computed")
    println("==========================================")
end

"""
    function compute_CoM!(obj)

Computes the integrals of motion of `obj.ψ` and saves them in respective fields of `obj`.
`obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_CoM`](@ref)
"""
function compute_CoM!(obj)
    bj
    if ~obj.spectrum_computed
        println("CoM calculation requested without a spectrum calculation.")
        compute_spectrum!(obj)
    end
    println("==========================================")
    println("Computing constants of motions")
    # Compute the energies
    # Do I need to find a better way of doing these integrals?
    # We should have 
    # sim.norm = sum(abs.(sim.ψ).^2, dims=2)[:]*sim.box.dt/sim.box.T but dt/T = 1/Nt, thus
    obj.N = sum(abs2.(obj.ψ), dims=2)[:]/obj.box.Nₜ
    obj.PE = -0.5*sum(abs2.(obj.ψ).^2,dims=2)[:]./(obj.N*obj.box.Nₜ)
    obj.KE = 0.5*sum((obj.box.ω'.^2) .* (abs2.(obj.ψ̃)),dims=2)[:]./obj.N
    obj.P = -imag.(sum(im * (obj.box.ω') .* (abs2.(obj.ψ̃)),dims=2)[:]./obj.N)
    obj.E = obj.KE + obj.PE
    obj.dE = obj.E .- obj.E[1]jjkkk
    println("Energy computed.")
    println("==========================================")
end

"""
    function ψ₀_periodic(coeff, box::SimBox, params::SimParamseters; phase=0)

Computes an initial wavefunction for the `SimBox` `box`, with fundamental frequency `sim.Ω`
and coefficients ``A_1...n`` = `coeff` and an overall phase `exp(i phase t)`, i.e. of the form:

**TODO**: insert latex form here

See also: [`init_sim`](@ref)
"""
function ψ₀_periodic(coeff::Array, box::Box, Ω; phase=0)
    println("==========================================")
    println("Initializing periodic ψ₀")
    for (n, An) in enumerate(coeff)
        if abs(An) >= 1
            println("The absolute value of the coefficient A($(n+1)) = $(An) is greater than 1. psi_0 needs to be normalizable.")
        else
            println("A($(n)) = $(An)")
        end #if
    end #for

    println("Computing A₀ to preserve normalization.")
    A0 = sqrt(1 - 2 * sum([abs(An) .^ 2 for An in coeff]))
    println("Computed A₀ = $A0")
    #A0 = 1
    if phase != 0
        str = "ψ₀ = exp(i $phase t) ($A0 + "
    else
        str = "ψ₀ = $A0 + "
    end

    ψ₀ = A0 * ones(box.Nₜ)
    for (n, An) in enumerate(coeff)
        ψ₀ += 2 * An * cos.(n * Ω * box.t)
        str = string(str, "2 × $An × cos($n × $(Ω) t)")
    end #for
    # Multiply by the overall phase
    ψ₀ = exp.(im * phase * box.t) .* ψ₀

    if phase != 0
        str = string(str, ")")
    end

    println(str)
    println("==========================================")

    return ψ₀
end #psi0_periodic