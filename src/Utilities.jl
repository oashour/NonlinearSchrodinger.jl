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
        # Why does the stuff from NumericalIntegration.jl not work?
        #obj.N .= integrate(obj.box.t, ψ², SimpsonEvenFast())./obj.T
        #obj.PE .= -0.5*integrate(obj.box.t, (ψ²).^2, SimpsonEvenFast())./obj.T
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
    println("Algorithm: $(sim.T̂)")
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

    return ψ₀, A0
end #psi0_periodic

function ψ₀_DT(λ, tₛ, xₛ, X₀, box)
    xᵣ = X₀=>X₀+1e-5
    T = abs(box.t[1]*2)
    Nₜ = box.Nₜ
    box_dt = Box(xᵣ, T, Nₓ=1, Nₜ = box.Nₜ)

    calc = Calc(λ, tₛ, xₛ, "exp", box_dt)
    solve!(calc)
    return calc.ψ[:, 1]

end

function λ_maximal(λ₁, N; m = 0)
    ν₁ = imag(λ₁)
    #ν_min = sqrt(1 - 1/N^2)
    mp = 1 - m
    C_N = sqrt((N^2-1)*mp*(N^2 - mp))
    H_N = N^2*(mp + 1) - 2*mp
    ν_min = sqrt(2*C_N + H_N)/(2*N)
    if ν₁ <= ν_min
        throw(ArgumentError("λ = $λ₁ not big enough for N = $N, need at least λ = $ν_min im"))
    end
    n = (1:N)
    #λ = sqrt.(n.^2 .* (ν₁^2 - 1) .+ 1)*im
    G_n = (m^2 .* n.^2) .+ 8*(m-2).*(n.^2 .- 1)*ν₁.^2 .+ 16 .* n.^2  .* ν₁^4
    λ = sqrt.(G_n .+ sqrt.(G_n.^2 .- 64*m^2*ν₁^4))./(4*sqrt(2)*ν₁)*im

    return λ
end