using NLSS
using Test

@testset "NLSS.jl" begin
    #TODO: Implement tests
    m = 0.2
    λ = λ_given_m(m, q=6)
    @test λ ≈ 0.9339432380535286im
    λ = λ_maximal(λ, 3, m=m)
    @test λ ≈ [0.9339432380535284im, 0.8929328208629054im, 0.819955288011106im]
    @test_throws ArgumentError λ_maximal(0.6im, 10, m=m)
end
