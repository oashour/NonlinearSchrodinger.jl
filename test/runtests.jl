using NLSS
using Test

@testset "NLSS.jl" begin
    params = NLSS.init_sim_params(a = 0.48)
    xᵣ = 0=>50
    box = NLSS.init_sim_box(xᵣ, params, Nₜ = 400, dx=1e-3, n_periods = 2)
    coeff = [1e-4]
    ψ₀ = NLSS.ψ₀_periodic(coeff, box, params)

    sim = NLSS.init_sim(box, params, ψ₀, algorithm = "4S", αₚ = 0)

    @time NLSS.solve!(sim)
    NLSS.plot_ψ(sim, mode="density")

    NLSS.compute_spectrum!(sim)
    NLSS.plot_ψ̃(sim, mode = "lines")

    NLSS.compute_energy!(sim)
    NLSS.plot_energy(sim)
end
