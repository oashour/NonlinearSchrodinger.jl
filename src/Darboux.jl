solve!(calc::Calc)
    ψₜ = Array{typeof{}}(undef, box.Nₜ, box.Nₓ)
end
function calc_rs(calc, n, p, ψₜ)

    # Base case
    if n == 1
        if calc.seed == "0"
            rf = exp(+im*(λ[p]*(x-xₛ[p]) + im*λ[p]^2*(t-tₛ[p]) - π/4))
            sf = exp(-im*(λ[p]*(x-xₛ[p]) + im*λ[p]^2*(t-tₛ[p]) - π/4))
        end
        if  p == 1
            if calc.seed == "0"
                ψ₀ = zeros(calc.box.Nₜ)
            end
            ψ[n] .= ψ₀ + (2*(λ[n]' - λ[n]).*sf.*conj.(rf))./(abs.(rf).^2 + abs.(sf).^2)
        end
    else
        # recursion
    end
end