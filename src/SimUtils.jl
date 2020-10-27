export Sim
export params, print
export ψ₀_periodic

###########################################################################
# Utility
###########################################################################
function params(; kwargs...)
    if length(kwargs) != 1
        throw(ArgumentError("You have either specified too few or too many parameters. You must specify one and only one of the following options: λ, Ω, T, a."))
    end
    param = Dict(kwargs)
    if :a in keys(param)
        λ = im * sqrt(2 * param[:a])
        T = π/sqrt(1 - imag(λ)^2)
        Ω = 2π/T
        @info "Passed a = $(param[:a]), computed λ = $λ, T = $T and Ω = $Ω"
    elseif :λ  in keys(param)
        λ = param[:λ]
        T = π/sqrt(1 - imag(λ)^2)
        Ω = 2π/T
        @info "Passed λ=$λ, computed T = $T and Ω = $Ω"
    elseif :Ω in keys(param)
        λ = im * sqrt((1 - (param[:Ω] / 2)^2))
        T = π/sqrt(1 - imag(λ)^2)
        Ω = param[:Ω]
        @info "Passed Ω=$Ω, computed λ = $λ and T = $T"
    elseif :T in keys(param)
        λ = im * sqrt((1 - ((2*π/param[:T]) / 2)^2))
        T = param[:T]
        Ω = 2π/T
        @info "Passed T = $T, computed λ = $λ and Ω = $Ω"
    end

    return λ, T, Ω
end #compute_parameters

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
    if sim.α == 0
        println("Equation: Cubic NLSE")
    elseif sim.α > 0
        println("Equation: Hirota Equation with α = $(sim.α)")
    else
        println("Unknown equation with α = $(sim.α)")
    end
    println("Algorithm of order $(sim.x_order) in x and $(sim.t_order) in t")
    println("------------------------------------------")
end

###########################################################################
# ψ₀
###########################################################################
"""
    function ψ₀_periodic(coeff, box::SimBox, params::SimParamseters; phase=0)

Computes an initial wavefunction for the `SimBox` `box`, with fundamental frequency `sim.Ω`
and coefficients ``A_1...n`` = `coeff` and an overall phase `exp(i phase t)`, i.e. of the form:

**TODO**: insert latex form here

See also: [`init_sim`](@ref)
"""
function ψ₀_periodic(coeff::Array, box::Box, Ω; phase=0)
    @info "Initializing periodic ψ₀"
    for (n, An) in enumerate(coeff)
        if abs(An) >= 1
            @error "The absolute value of the coefficient A($(n+1)) = $(An) is greater than 1. psi_0 needs to be normalizable."
        end #if
    end #for

    @info "Computing A₀ to preserve normalization."
    A0 = sqrt(1 - 2 * sum([abs(An) .^ 2 for An in coeff]))
    @info "Computed A₀ = $A0"
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

    @info str
    println("==========================================")

    return ψ₀
end #psi0_periodic