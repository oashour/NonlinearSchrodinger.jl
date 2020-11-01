function solve!(calc::Calc)
    N = length(calc.λ)
    ψₜ = Array{Complex{Float64}}(undef, calc.box.Nₜ, calc.box.Nₓ, N+1)
    if calc.seed == "0"
        ψₜ[:,:,1] .= zeros(calc.box.Nₜ, calc.box.Nₓ)
    elseif calc.seed == "exp"
        ψₜ[:,:,1] .= exp.(im*calc.box.x').*ones(calc.box.Nₜ, calc.box.Nₓ)
    elseif calc.seed == "dn"
        ψₜ[:,:,1] .= exp.(im*calc.box.x'.*(1-calc.m/2)).*dn.(calc.box.t,calc.m)
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
        elseif c.seed == "exp"
            t = c.box.t .- c.tₛ[p]
            x = c.box.x' .- c.xₛ[p]
            A = +c.χ[p] .+ 0.5.*(c.Ω[p].*t .+ c.Ω[p]*c.λ[p].*x) .- π/4
            B = -c.χ[p] .+ 0.5.*(c.Ω[p].*t .+ c.Ω[p]*c.λ[p].*x) .- π/4
            rf = 2im.*exp.(-im.*c.box.x'./2) .* sin.(A) # why not the shifted x?
            sf = 2  .*exp.(+im.*c.box.x'./2) .* cos.(B) # why not the shifted x?
        elseif c.seed == "dn"
            rf = Array{Complex{Float64}}(undef, c.box.Nₜ, c.box.Nₓ)
            sf = Array{Complex{Float64}}(undef, c.box.Nₜ, c.box.Nₓ)
            #t = c.box.t .- c.tₛ[p]
            x = c.box.x' .- c.xₛ[p]
            function dn_ab!(du,u,p,t)
                du[1] = +im*p[1].*u[1] .+ im.*u[2].*dn.(t,real(p[2]))
                du[2] = -im*p[1].*u[2] .+ im.*u[1].*dn.(t,real(p[2]))
            end

            function a_dn_t0(x)
                A = +c.χ[p] .+ 0.5.*c.Ω[p]*c.λ[p].*(x - c.xₛ[p]) .- π/4
                2*im*sin(A)
            end
            function b_dn_t0(x)
                B = -c.χ[p] .+ 0.5.*c.Ω[p]*c.λ[p].*(x - c.xₛ[p]) .- π/4
                2*cos(B)
            end
            tspan = (0.0,abs(minimum(c.box.t)))
            for i in eachindex(c.box.x)
                u0 = [a_dn_t0(c.box.x[i]);b_dn_t0(c.box.x[i])]
                prob = ODEProblem(dn_ab!, u0, tspan, [c.λ[p], c.m])
                sol = solve(prob, Tsit5(), saveat=abs.(c.box.t[c.box.Nₜ÷2+1:-1:1]))
                rtemp = sol[1, :].*exp.(+im*(c.box.x[i] - c.xₛ[p])/4*(c.m-2))
                stemp = sol[2, :].*exp.(-im*(c.box.x[i] - c.xₛ[p])/4*(c.m-2))
                rf[c.box.Nₜ÷2+1:-1:1, i] = rtemp
                sf[c.box.Nₜ÷2+1:-1:1, i] = stemp
                rf[c.box.Nₜ÷2+2:end, i] = rtemp[2:end-1] 
                sf[c.box.Nₜ÷2+2:end, i] = stemp[2:end-1] 
            end
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