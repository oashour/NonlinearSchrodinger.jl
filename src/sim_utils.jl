export print, ψ₀_periodic, compute_energy!, compute_spectrum!

function print(sim::Simulation)
    println("Box Properties:")
    println("------------------------------------------")
    println("dx = $(sim.box.dx) (Nₓ = $(sim.box.Nₓ))")
    println("Nₜ = $(sim.box.dx) (dt = $(sim.box.dt))")
    println("Parameters:")
    println("------------------------------------------")
    println("a = $(sim.params.a)")
    println("λ = $(sim.params.λ)")
    println("Ω = $(sim.params.Ω)")
    println("T = $(sim.params.T)")
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

function compute_spectrum!(sim::Simulation)
    println("==========================================")
    println("Computing spectrum")
    if sim.solved
        sim.ψ̃ = fftshift(fft(sim.ψ, 2), 2)/sim.box.Nₜ
        sim.spectrum_computed = true
    else
        throw(ArgumentError("Trying to compute spectrum of an unsolved simulation. Please solve the model first."))
    end
    println("Spectrum computed")
    println("==========================================")
end

function compute_CoM!(sim::Simulation)
    if ~sim.spectrum_computed
        println("CoM calculation requested without a spectrum calculation.")
        compute_spectrum!(sim)
    end
    println("==========================================")
    println("Computing constants of motions")
    # Compute the energies
    # Do I need to find a better way of doing these integrals?
    # We should have 
    # sim.norm = sum(abs.(sim.ψ).^2, dims=2)[:]*sim.box.dt/sim.box.T but dt/T = 1/Nt, thus
    sim.N = sum(abs2.(sim.ψ), dims=2)[:]/sim.box.Nₜ
    sim.PE = -0.5*sum(abs2.(sim.ψ).^2,dims=2)[:]./(sim.N*sim.box.Nₜ)
    sim.KE = 0.5*sum((sim.box.ω'.^2) .* (abs2.(sim.ψ̃)),dims=2)[:]./sim.N
    sim.P = -imag.(sum(im * (sim.box.ω') .* (abs2.(sim.ψ̃)),dims=2)[:]./sim.N)
    sim.E = sim.KE + sim.PE
    sim.dE = sim.E .- sim.E[1]
    println("Energy computed.")
    println("==========================================")
end

function ψ₀_periodic(coeff, box::SimulationBox, params::SimulationParameters; phase=0)
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
        ψ₀ += 2 * An * cos.(n * params.Ω * box.t)
        str = string(str, "2 × $An × cos($n × $(params.Ω) t)")
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