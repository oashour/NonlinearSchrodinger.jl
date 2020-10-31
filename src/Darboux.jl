function solve!(calc::Calc)
    N = length(calc.λ)
    ψₜ = Array{Complex{Float64}}(undef, calc.box.Nₜ, calc.box.Nₓ, N+1)
    if calc.seed == "0"
        ψₜ[:,:,1] .= zeros(calc.box.Nₜ, calc.box.Nₓ)
    end
    calc_rs(calc, N, 1, ψₜ)
    calc.ψ .= ψₜ[:,:,N+1]
    calc.ψ̃ .= fftshift(fft(calc.ψ, 1),1)./calc.box.Nₜ
end
function calc_rs(c::Calc, n, p, ψₜ)
    @info "Calculating Lax pair generating functions rₙₚ(x,t) and sₙₚ(x,t) for (n,p) = ($n,$p)"

    # Base case
    if n == 1
        if c.seed == "0"
            t = c.box.t .- c.tₛ[p]
            x = c.box.x' .- c.xₛ[p]
            rf = exp.(+im.*(c.λ[p] .* t .+ c.λ[p]^2 .* x .- π/4))
            sf = exp.(-im.*(c.λ[p] .* t .+ c.λ[p]^2 .* x .- π/4))
        end
        if  p == 1
            ψₜ[:,:,n+1] .= ψₜ[:,:,1] .+ (2*(c.λ[n]' - c.λ[n]).*sf.*conj.(rf)) ./
                              (abs2.(rf) + abs2.(sf))
        end
    else
        r1, s1 = calc_rs(c, n-1, 1, ψₜ)
        r2, s2 = calc_rs(c, n-1, p+1, ψₜ)

        rf = ((c.λ[n-1]'  - c.λ[n-1] ) .* conj.(s1) .* r1 .* s2 +
              (c.λ[p+n-1] - c.λ[n-1] ) .* abs2.(r1) .* r2       +
              (c.λ[p+n-1] - c.λ[n-1]') .* abs2.(s1) .* r2        ) ./
              (abs.(r1).^2 + abs.(s1).^2);
        sf = ((c.λ[n-1]'  - c.λ[n-1] ) .*       s1   .*conj(r1) .* r2 +
              (c.λ[p+n-1] - c.λ[n-1] ) .*  abs2.(s1) .*     s2        +
              (c.λ[p+n-1] - c.λ[n-1]') .*  abs2.(r1) .*     s2         ) ./
              (abs.(r1).^2 + abs.(s1).^2);
        if p == 1
            ψₜ[:,:,n+1] .= ψₜ[:,:,n] + (2*(c.λ[n]' - c.λ[n]).*sf.*conj.(rf)) ./
                                      (abs2.(rf) + abs2.(sf))
        end
    end

    return rf, sf
end