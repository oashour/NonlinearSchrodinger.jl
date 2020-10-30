function solve!(calc::Calc)
    N = length(calc.λ)
    ψₜ = Array{Complex{Float64}}(undef, calc.box.Nₜ, calc.box.Nₓ, N)
    calc_rs(calc, N, 1, ψₜ)
    calc.ψ .= ψₜ[:,:,N]
end
function calc_rs(calc, n, p, ψₜ)

    # Base case
    if n == 1
        if calc.seed == "0"
            rf = exp.(+im.*(calc.λ[p]  .*(calc.box.t  .- calc.tₛ[p]) .+ 
                           calc.λ[p]^2 .*(calc.box.x' .- calc.xₛ[p]) .- π/4))
            sf = exp.(-im.*(calc.λ[p]  .*(calc.box.t  .- calc.tₛ[p]) .+ 
                           calc.λ[p]^2 .*(calc.box.x' .- calc.xₛ[p]) .- π/4))
        end
        if  p == 1
            if calc.seed == "0"
                ψ₀ = zeros(calc.box.Nₜ, calc.box.Nₓ)
            end
            ψₜ[:,:,n] .= ψ₀ .+ (2*(calc.λ[n]' - calc.λ[n]).*sf.*conj.(rf)) ./
                              (abs2.(rf) + abs2.(sf))
        end
    else
        r1, s1 = calc_rs(calc, n-1, 1, ψₜ)
        r2, s2 = calc_rs(calc, n-1, p+1, ψₜ)

        rf = ((calc.λ[n-1]'  - calc.λ[n-1] ) .* conj.(s1) .* r1 .* s2 +
              (calc.λ[p+n-1] - calc.λ[n-1] ) .* abs2.(r1) .* r2       +
              (calc.λ[p+n-1] - calc.λ[n-1]') .* abs2.(s1) .* r2        ) ./
              (abs.(r1).^2 + abs.(s1).^2);
        sf = ((calc.λ[n-1]'  - calc.λ[n-1] ) .*       s1   .*conj(r1) .* r2 +
              (calc.λ[p+n-1] - calc.λ[n-1] ) .*  abs2.(s1) .*     s2        +
              (calc.λ[p+n-1] - calc.λ[n-1]') .*  abs2.(r1) .*     s2         ) ./
              (abs.(r1).^2 + abs.(s1).^2);
        if p == 1
            ψₜ[:,:,n] .= ψₜ[:,:,n-1] + (2*(calc.λ[n]' - calc.λ[n]).*sf.*conj.(rf)) ./
                                      (abs2.(rf) + abs2.(sf))
        end
    end

    return rf, sf
end