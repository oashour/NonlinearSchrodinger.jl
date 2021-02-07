using NonlinearSchrodinger
using Test

@testset "NonlinearSchrodinger.jl" begin
    # Tests on maximal intensity generating functions
    m = 0.2
    λ = λ_given_m(m, q=6)
    @test λ ≈ 0.9339432380535286im
    λ = λ_maximal(λ, 3, m=m)
    @test λ ≈ [0.9339432380535284im, 0.8929328208629054im, 0.819955288011106im]
    @test_throws ArgumentError λ_maximal(0.6im, 10, m=m)

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
                T4A_N!, T6A_N!, T8A_N!, T6A_OP!, T8A_OP!,
                T2B!, T4B_TJ!, T6B_TJ!, T8B_TJ!, T4B_SF!, T4B_SF!, T6B_SF!, T8B_SF!,
                T4B_N!, T6B_N!, T8B_N!, T6B_OP!, T8B_OP!)
            sim = Sim(λ, box, ψ₀, algo)
            @info "Testing Algorithm" algo
            solve!(sim)
            @test abs.(sim.ψ[:, end]) ≈ abs.(ψ₀) atol=1e-5
        end
    end
end
