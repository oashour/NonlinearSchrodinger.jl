using NonlinearSchrodinger
using RecipesBase
using Test

@testset "NonlinearSchrodinger.jl" begin
    # Tests on maximal intensity generating functions
    m = 0.2
    λ = λ_given_m(m, q=6)
    @test λ ≈ 0.9339432380535286im
    λ = λ_maximal(λ, 3, m=m)
    @test λ ≈ [0.9339432380535284im, 0.8929328208629054im, 0.819955288011106im]
    @test_throws ArgumentError λ_maximal(0.6im, 10, m=m)

    @testset "Initial Conditions" begin
        xᵣ = -10=>10
        λ = 0.75im
        λ, T, Ω = params(λ = λ)
        box = Box(xᵣ, T, Nₓ=1, Nₜ = 100, n_periods = 3)
        coeff = [(2.7 + 4.6im)*1e-2]
        ψ₀, A₀ = ψ₀_periodic(coeff, box, Ω)
        @test A₀ ≈ 0.9971509414326399
    end

    @testset "params" begin
        λ, T, Ω = params(a=3/8)
        @test imag(λ) ≈ sqrt(3)/2
        @test T ≈ 2*π
        @test Ω ≈ 1
        λ, T, Ω = params(λ = sqrt(3)/2*im)
        @test T ≈ 2*π
        @test Ω ≈ 1
        λ, T, Ω = params(T = 2*pi) 
        @test imag(λ) ≈ sqrt(3)/2
        @test Ω ≈ 1
        λ, T, Ω = params(Ω = 1) 
        @test imag(λ) ≈ sqrt(3)/2
        @test T ≈ 2*π
    end

    @testset "Darboux Transformations" begin
        @testset "Soliton 1" begin
            T = 20
            xᵣ = -20=>20
            seed = "0"
            box = Box(xᵣ, T, Nₓ=1000, Nₜ = 1024, n_periods = 1)
            λ = [0.85im]
            xₛ = [0.0]
            tₛ = [0.0]

            calc = Calc(λ, tₛ, xₛ, seed, box) 
            solve!(calc)
            λ = 0.85im; 
            ψ₀ = 2*imag(λ)./cosh.(2*imag(λ).*box.t)
            @test abs.(calc.ψ[:, end]) ≈ abs.(ψ₀)
        end
        @testset "Soliton 3" begin
            T = 20
            xᵣ = -20=>20
            seed = "0"
            box = Box(xᵣ, T, Nₓ=1000, Nₜ = 1024, n_periods = 1)
            λ = [-0.9 + 0.8im, 0.85im, 0.9 + 0.9im]
            xₛ = [0.0, 0.0, 0.0]
            tₛ = [0.0, 0.0, 0.0]

            calc = Calc(λ, tₛ, xₛ, seed, box) 
            solve!(calc)
            peak = PHF(calc)
            @test peak ≈ 5.1

            λ = 0.85im; 
            ψ₀ = 2*imag(λ)./cosh.(2*imag(λ).*box.t)
            # Test around edges and exact peak to avoid weird behavior near peak from DT iteration
            @test abs.(calc.ψ[1:250, end]) ≈ abs.(ψ₀[1:250]) atol=1e-3
            @test abs.(calc.ψ[750:end, end]) ≈ abs.(ψ₀[750:end]) atol=1e-3
            @test abs.(calc.ψ[513, end]) ≈ abs.(ψ₀[513]) atol=1e-2
        end
        @testset "AB 1" begin
            T = 20
            xᵣ = 0=>1e-3
            seed = "exp"
            box = Box(xᵣ, T, Nₓ=1, Nₜ = 1024, n_periods = 1)
            λ = [0.9im]
            xₛ = [0.0]
            tₛ = [0.0]

            calc = Calc(λ, tₛ, xₛ, seed, box) 
            solve!(calc)
            peak = PHF(calc)
            @test peak ≈ 2.8

            λ = 0.9im; 
            a = imag(λ)^2/2
            Ω = 2*sqrt(1 - 2*a)
            ψ₀ = ((1 - 4*a) .+ (sqrt(2*a).*cos.(Ω*calc.box.t)))./
                (sqrt(2*a).*cos.(Ω*calc.box.t) .- 1)

            @test abs.(calc.ψ[:, 1]) ≈ abs.(ψ₀)
        end
        @testset "CN 1" begin
            T = 20
            xᵣ = -10=>10
            box = Box(xᵣ, T, Nₓ=1, Nₜ = 1024, n_periods = 1)
            λ = [0.7im]
            xₛ = [0.0]
            tₛ = [0.0]

            seed = "cn"
            calc = Calc(λ, tₛ, xₛ, seed, box, m = 0.5) 
            peak = PHF(calc)
            @test peak ≈ 2.1071067811865474

            solve!(calc)

            # Test the solution is mirrored across the x=0 axis
            @test abs.(calc.ψ[:,2]) ≈ abs.(calc.ψ[:,1])
        end
        @testset "DN 1" begin
            T = 20
            xᵣ = -10=>10
            box = Box(xᵣ, T, Nₓ=1, Nₜ = 1024, n_periods = 1)
            λ = [0.6im]
            xₛ = [0.0]
            tₛ = [0.0]

            seed = "dn"
            calc = Calc(λ, tₛ, xₛ, seed, box, m = 0.5) 
            peak = PHF(calc)
            @test peak ≈ 2.2

            solve!(calc)

            # Test the solution is mirrored across the x=0 axis
            @test abs.(calc.ψ[:,2]) ≈ abs.(calc.ψ[:,1])
        end
    end

    # Simulation Tests using Solitons (easiest way to check validity)
    @testset "Soliton Simulations" begin
        λ = 0.75im
        T = 20
        xᵣ = 0=>10
        box = Box(xᵣ, T, dx=1e-3, Nₜ = 256, n_periods = 1)
        ψ₀ = Array{Complex{Float64}}(undef, box.Nₜ)
        ψ₀ .= 2*imag(λ)./cosh.(2*imag(λ).*box.t)

        for algo ∈ (T1A!, T1B!)
            sim = Sim(λ, box, ψ₀, algo)
            @info "Testing Algorithm" algo
            solve!(sim)
            @test abs.(sim.ψ[:, end]) ≈ abs.(ψ₀) atol=1e-2
        end
        for algo ∈ (T2A!, T4A_TJ!, T6A_TJ!, T8A_TJ!, T4A_SF!, T4A_SF!, T6A_SF!, T8A_SF!, 
                T4A_CMP!, T6A_CMP!, T8A_CMP!, T6A_Ss14!, T8A_Ss15!,
                T2B!, T4B_TJ!, T6B_TJ!, T8B_TJ!, T4B_SF!, T4B_SF!, T6B_SF!, T8B_SF!,
                T4B_CMP!, T6B_CMP!, T8B_CMP!, T6B_Ss14!, T8B_Ss15!)
            sim = Sim(λ, box, ψ₀, algo)
            @info "Testing Algorithm" algo
            solve!(sim)
            @test abs.(sim.ψ[:, end]) ≈ abs.(ψ₀) atol=1e-5
        end
    end
    @testset "Soliton Simulations Hirota" begin
        λ = 0.75im
        T = 20
        xᵣ = 0=>5
        box = Box(xᵣ, T, dx=1e-3, Nₜ = 256, n_periods = 1)
        ψ₀ = Array{Complex{Float64}}(undef, box.Nₜ)
        ψf = Array{Complex{Float64}}(undef, box.Nₜ)
        ψ₀ .= 2*imag(λ)./cosh.(2*imag(λ).*box.t)
        xf = box.x[end]
        α = 0.1
        v = 4*α*imag(λ)^2
        ψf .= 2*imag(λ)./cosh.(2*imag(λ).*(box.t .+ v*xf))

        for algo ∈ (T1A_H!, T2A_H!)
            sim = Sim(λ, box, ψ₀, algo, α = α)
            @info "Testing Algorithm" algo
            solve!(sim)
            # Can't test around the peak due to a minor shift between them
            @test abs.(sim.ψ[1:50, end]) ≈ abs.(ψ₀[1:50]) atol=5e-2
            @test abs.(sim.ψ[175:end, end]) ≈ abs.(ψ₀[175:end]) atol=5e-2
        end
    end
    @testset "Plot Recipes" begin
        λ = 0.75im
        T = 20
        xᵣ = 0=>1
        box = Box(xᵣ, T, dx=1e-3, Nₜ = 256, n_periods = 1)
        ψ₀ = Array{Complex{Float64}}(undef, box.Nₜ)
        ψ₀ .= 2*imag(λ)./cosh.(2*imag(λ).*box.t)
        sim = Sim(λ, box, ψ₀, T1A!)
        solve!(sim)
        compute_IoM!(sim)
        xₛ = Int(ceil(sim.box.Nₓ/500))
        tₛ = Int(ceil(sim.box.Nₜ/500))
        x_ax = sim.box.t[1:tₛ:end]
        y_ax =  sim.box.x[1:xₛ:end] 
        ω_ax = sim.box.ω[1:tₛ:end] 
        ψ = sim.ψ[1:tₛ:end, 1:xₛ:end]
        ψ̃ = sim.ψ̃[1:tₛ:end, 1:xₛ:end]

        # Test plotting ψ
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), sim, :ψ)
        @test getfield(rec[1], 2)[1] == x_ax
        @test getfield(rec[1], 2)[2] == y_ax
        @test getfield(rec[1], 2)[3] == abs.(ψ')

        # Test plotting ψ̃
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), sim, :ψ̃)
        @test getfield(rec[1], 2)[1] == y_ax
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(:seriestype => :heatmap), sim, :ψ̃)
        @test getfield(rec[1], 2)[1] == ω_ax
        @test getfield(rec[1], 2)[2] == y_ax
        @test getfield(rec[1], 2)[3] == log.(abs.(ψ̃'))

        # Test plotting IoM
        rec = RecipesBase.apply_recipe(Dict{Symbol, Any}(), sim, :IoM)
        @test getfield(rec[1], 2)[1] == y_ax
        @test getfield(rec[2], 2)[1] == y_ax
        @test getfield(rec[3], 2)[1] == y_ax
        @test getfield(rec[4], 2)[1] == y_ax

        # Test Error
        @test_logs (:error, "Unknown Plotting Mode for NonlinearSchrodinger.jl Objects IoMM") RecipesBase.apply_recipe(Dict{Symbol, Any}(), sim, :IoMM)
         
    end
end
