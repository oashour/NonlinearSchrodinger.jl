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
    @info "Initializing simulation box with $n_periods period(s) and dx = $dx, Nₜ = $Nₜ."
    T = n_periods * T
    println("Longitudinal range is [$(xᵣ.first), $(xᵣ.second)], transverse range is [$(-T/2), $(T/2))")
    dt = T / Nₜ
    t = dt * collect((-Nₜ/2:Nₜ/2-1))

    x = collect(xᵣ.first:dx:xᵣ.second)
    Nₓ = length(x)
    ω = 2π/T * collect((-Nₜ/2:Nₜ/2-1))

    box = Box(t, ω, x, Nₜ, Nₓ, dt, dx, n_periods)

    @info "Done computing t, x, ω"

    return box
end

"""
    function compute_spectrum!(obj)

Computes the normalized spectrum of `obj.ψ` with the center frequency shifted to the center
and saves it in `obj.ψ̃`. `obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_ψ̃`](@ref)
"""
function compute_spectrum!(obj)
    @info "Computing spectrum"

    FFTW.set_num_threads(4)
    F̂ = plan_fft(obj.ψ, 1) # 26 allocs
    obj.ψ̃ .= fftshift(F̂*obj.ψ, 1)/obj.box.Nₜ
    FFTW.set_num_threads(1) # set it back to 1, not useful in smaller FFTs

    @info "Spectrum computed"
end

"""
    function compute_IoM!(obj)

Computes the integrals of motion of `obj.ψ` and saves them in respective fields of `obj`.
`obj` can be a `Sim` or `Calc` object

See also: [`NLSS.Plotter.plot_CoM`](@ref)
"""
function compute_IoM!(obj; normalize=false)
    @info "Computing integrals of motions"
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

    @info "Integrals of motion computed."
end

function save(obj, filename)
    filename = string(filename, ".jld")
    jldopen(filename, "w") do file
        write(file, "result", obj)
    end
    @info "Saved to file $filename"
    
    return nothing
end

function load(filename)
    filename = string(filename, ".jld")
    obj = jldopen(filename, "r") do file
        read(file, "result")
    end
    @info "Loaded from file $filename"

    return obj
end

"""
    function params(; kwargs...)

Computes parameters
"""
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