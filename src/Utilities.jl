"""
    function compute_IoM!(obj)

Computes the integrals of motion of `obj.ψ` and saves them in respective fields of `obj`.
`obj` can be a `::Sim` or `::Calc` object. The fields are:

```
obj.KE  # Array containing the kinetic energy K(x)
obj.PE  # Array containing the potential energy V(x)
obj.E   # Array containing the energy H(x)
obj.N   # Array containing the norm N(x)
obj.P   # Array containing the momentum P(x)
obj.dE  # Array containing the energy error δE(x)
obj.dN  # Array containing the norm error δN(x)
obj.dP  # Array containing the momentum error δP(x)
```

Result can be plotted using `plot(obj, :IoM)`
"""
function compute_IoM!(obj)
    @info "Computing integrals of motions"
    ψ² = abs2.(obj.ψ)
    ψ̃² = abs2.(obj.ψ̃)
    # The norm is preserved extremely well so if the W.F. is normalized we do not worry about dividing the IoM by the norm to improve performance
    #if normalize
        #obj.N .= sum(ψ², dims=1)[:]./obj.box.Nₜ
        #obj.PE .= -0.5*sum(ψ².^2,dims=1)[:]./(obj.N*obj.box.Nₜ)
        #obj.KE .= 0.5*sum((obj.box.ω.^2 .* ψ̃²),dims=1)[:]./obj.N
        #obj.P .= 2*imag.(sum(im * obj.box.ω.* ψ̃²,dims=1)[:]./obj.N)
        #obj.P .= 2*sum(obj.box.ω.* ψ̃²,dims=1)[:]./obj.N
        #obj.E .= obj.KE .+ obj.PE
    #else
    # Why does the stuff from NumericalIntegration.jl not work?
    #obj.N .= integrate(obj.box.t, ψ², SimpsonEvenFast())./obj.T
    #obj.PE .= -0.5*integrate(obj.box.t, (ψ²).^2, SimpsonEvenFast())./obj.T
    obj.N .= sum(ψ², dims=1)[:]./obj.box.Nₜ
    obj.PE .= -0.5*sum(ψ².^2,dims=1)[:]./obj.box.Nₜ
    obj.KE .= 0.5*sum((obj.box.ω.^2 .* ψ̃²),dims=1)[:]
    #obj.P .= 2*imag.(sum(im * obj.box.ω.* ψ̃²,dims=1)[:])
    obj.P .= 2*sum(obj.box.ω.* ψ̃²,dims=1)[:]
    obj.E .= obj.KE .+ obj.PE
    #end

    @info "Integrals of motion computed."
end

"""
    function params(; kwargs...)

Computes parameters `λ, T, Ω` given either of `λ, T, Ω, a` and possibly `m` for dnoidal seed. Used only for `exp` and `dn` seeds.
"""
function params(;m=0.0, kwargs...)
    if length(kwargs) != 1
        throw(ArgumentError("You have either specified too few or too many parameters. You must specify one and only one of the following options: λ, Ω, T, a."))
    end
    param = Dict(kwargs)
    if m ==0.0
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
    elseif m > 0 && m <= 1
        if :a in keys(param)
            @error "Passing a not yet supported in m != 0 mode. Please pass λ"
        elseif :λ  in keys(param)
            λ = param[:λ]
            Ω = real(2*sqrt(1 + (λ - m/4/λ)^2))
            T = 2π/Ω
            @info "Passed λ=$λ, computed T = $T and Ω = $Ω"
        elseif :Ω in keys(param)
            @error "Passing Ω not yet supported in m != 0 mode. Please pass λ"
        elseif :T in keys(param)
            @error "Passing T not yet supported in m != 0 mode. Please pass λ"
        end
    elseif m > 1 || m < 0
        @error "Wrong range for m. m must be between 0 and 1"
    end

    return λ, T, Ω
end #compute_parameters

###########################################################################
# ψ₀
###########################################################################
"""
    function ψ₀_periodic(coeff::Array, box::Box, Ω; phase=0)

Computes an initial wavefunction for the `box::Box`, with fundamental frequency `Ω`
and coefficients ``A_{1\\ldots n}`` = `coeff` and an overall phase ``e^{i \\phi t}`` where ``\\phi`` = `phase`. i.e. ``\\psi_0`` is of the form:

``
\\psi(x=0, t) = e^{i \\phi t} (A_0 + 2 \\sum_{1}^{n} A_m \\cos(m \\Omega t))
``

where:

``
A_0 = \\sqrt{1 - 2 \\sum_{m=1}^n |A_m|^2}
``
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

    return ψ₀, A0
end #psi0_periodic

"""
    ψ₀_DT(λ, tₛ, xₛ, X₀, box; seed="exp", f = Dict(:α=> 0.0, :γ => 0.0, :δ=>0.0))

Computes an initial condition using the DT characterized by `λ, xₛ, tₛ` and `seed` at x = `X₀` for the extended equation characterized by a dictionary of parameters `f` (currently only `α`) is supported. `box::Box` is the simulation box.

See also: [`Calc`](@ref), [`Sim`](@ref)
"""
function ψ₀_DT(λ, tₛ, xₛ, X₀, box; seed="exp", f = Dict(:α=> 0.0, :γ => 0.0, :δ=>0.0))
    xᵣ = X₀=>X₀+1e-5
    T = abs(box.t[1]*2)
    Nₜ = box.Nₜ
    box_dt = Box(xᵣ, T, Nₓ=1, Nₜ = box.Nₜ)

    calc = Calc(λ, tₛ, xₛ, seed, box_dt, f = f)
    solve!(calc)

    return calc.ψ[:, 1]
end

"""
    λ_maximal(λ₁, N; m = 0)

Computed a maximal intensity set of `λ` of order `N` given `λ₁` and possibly `m` for dnoidal bakcground.

See also: [`λ_given_m`](@ref)
"""
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
end

"""
    λ_given_m(m; q = 2)

Computed a `λ` that is matched to the dnoidal background given by `m` with an integer `q`. See paper for more details.

See also: [`λ_maximal`](@ref)
"""
function λ_given_m(m; q = 2)
    F = π/(2*q*Elliptic.K(m))
    λ = 0.5 * sqrt(2 - 2*F^2 - m + 2*sqrt((F^2-1)*(F^2-1+m)))*im
end

"""
    λ_given_f(f, ν)

Computes `λ = v + i ν` such that `v` results in breather to soliton conversion for the extended NLS characterized by the dictionary of parameters `f`.
"""
function λ_given_f(f, ν)
    α = f[:α]
    γ = f[:γ]
    δ = f[:δ]
    @info "Computing real part of eigenvalue for breather to soliton conversion"
    @info "Got α = $α, γ = $γ, δ = $δ"
    poly = reverse([64*δ, -24*γ, -8*(α+2*δ+8*δ*ν^2), γ*(4 + 8*ν^2)+1])
    r = roots(poly)

    if imag(r[1]) == 0
        @info "Computed v = $(real(r[1]))"
        return real(r[1]) + ν*im
    elseif imag(r[2]) == 0
        @info "Computed v = $(real(r[2]))"
        return real(r[2]) + ν*im
    elseif imag(r[3]) == 0
        @info "Computed v = $(real(r[3]))"
        return real(r[3]) + ν*im
    end
end

"""
    PHF(calc::Calc)

Computed the peak height formula for any `calc::Calc` objecti irrespective of seed.
"""
function PHF(calc::Calc)
    s = 2*sum(imag.(calc.λ))
    if calc.seed == "exp" || calc.seed == "dn"
        ψ₀₀ = 1
    elseif calc.seed == "cn"
        ψ₀₀ = sqrt(calc.m)
    elseif calc.seed == "0"
        ψ₀₀ = 0
    end
    peak = ψ₀₀ + s
end

