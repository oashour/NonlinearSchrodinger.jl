function solve!(calc::Calc)
    N = length(calc.λ)
    ψₜ = Array{Complex{Float64}}(undef, calc.box.Nₜ, calc.box.Nₓ, N+1)
    if calc.seed == "0"
        ψₜ[:,:,1] .= zeros(calc.box.Nₜ, calc.box.Nₓ)
    elseif calc.seed == "exp"
        ψₜ[:,:,1] .= exp.(im*calc.box.x').*ones(calc.box.Nₜ, calc.box.Nₓ)
    elseif calc.seed == "dn"
        ψₜ[:,:,1] .= exp.(im*calc.box.x'.*(1-calc.m/2)).*dn.(calc.box.t,calc.m)
    elseif calc.seed == "cn"
        ψₜ[:,:,1] .= exp.(im*calc.box.x'.*(calc.m - 1/2)).*sqrt(calc.m).*cn.(calc.box.t,calc.m)
    end
    calc_rs(calc, N, 1, ψₜ)
    calc.ψ .= ψₜ[:,:,N+1]
    calc.ψ̃ .= fftshift(fft(calc.ψ, 1),1)./calc.box.Nₜ
end

@memoize function calc_rs(c::Calc, n, p, ψₜ)

    @info "Calculating Lax pair generating functions rₙₚ(x,t) and sₙₚ(x,t) for (n,p) = ($n,$p)"
    # Base case
    if n == 1
        if c.seed == "0"
            # Set up x and t in the proper way
            t = c.box.t .- c.tₛ[p]
            x = c.box.x' .- c.xₛ[p]
            # Compute r and s directly
            rf = exp.(+im.*(c.λ[p] .* t .+ c.λ[p]^2 .* x .- π/4))
            sf = exp.(-im.*(c.λ[p] .* t .+ c.λ[p]^2 .* x .- π/4))
        elseif c.seed == "exp"
            # Set up x and t in the proper way
            t = c.box.t .- c.tₛ[p]
            x = c.box.x' .- c.xₛ[p]
            # Compute r and s directly
            A = +c.χ[p] .+ 0.5.*(c.Ω[p].*t .+ c.Ω[p]*c.λ[p].*x) .- π/4
            B = -c.χ[p] .+ 0.5.*(c.Ω[p].*t .+ c.Ω[p]*c.λ[p].*x) .- π/4
            rf = 2im.*exp.(-im.*c.box.x'./2) .* sin.(A)
            sf = 2  .*exp.(+im.*c.box.x'./2) .* cos.(B)
        elseif c.seed == "dn"
            # Allocate empty arrays to hold lax pair generating functions
            rf = Array{Complex{Float64}}(undef, c.box.Nₜ, c.box.Nₓ)
            sf = Array{Complex{Float64}}(undef, c.box.Nₜ, c.box.Nₓ)
            # Set up x in the proper way
            # We currently do not suppport t shifts for this seed
            x = c.box.x' .- c.xₛ[p]

            # Set up a funnction for the ODE
            function dn_ab!(du,u,p,t)
                du[1] = +im*p[1].*u[1] .+ im.*u[2].*dn.(t,real(p[2]))
                du[2] = -im*p[1].*u[2] .+ im.*u[1].*dn.(t,real(p[2]))
            end
            # And the initial condition
            function ab_dn_t0(x)
                A = +c.χ[p] .+ 0.5.*c.Ω[p]*c.λ[p].*(x - c.xₛ[p]) .- π/4
                B = -c.χ[p] .+ 0.5.*c.Ω[p]*c.λ[p].*(x - c.xₛ[p]) .- π/4
                return [2*im*sin(A); 2*cos(B)]
            end

            # Prepare the ODE intnegrator
            tspan = (0.0,abs(minimum(c.box.t)))
            u0 = ab_dn_t0(c.box.x[1])
            prob = ODEProblem(dn_ab!, u0, tspan, [c.λ[p], c.m])
            dt = c.box.t[2]-c.box.t[1]
            integrator = init(prob, Tsit5(), saveat=abs.(c.box.t[c.box.Nₜ÷2+1:-1:1]), dt = dt, adaptive=false)

            # Loop over x and solve the ODE for each initial condition
            # Should do this using the ensemble interface in the future
            for i in eachindex(c.box.x)
                # Set up the new initial condition
                u0 = ab_dn_t0(c.box.x[i])
                reinit!(integrator, u0)
                # Solve
                DiffEqBase.solve!(integrator)
                # Save results
                rf[c.box.Nₜ÷2+1:-1:1, i] .= integrator.sol[1, :].*exp.(+im*(c.box.x[i])./4*(c.m-2))
                sf[c.box.Nₜ÷2+1:-1:1, i] .= integrator.sol[2, :].*exp.(-im*(c.box.x[i])./4*(c.m-2))
                # Mirror x -> x
                rf[c.box.Nₜ÷2+2:end, i] .= rf[c.box.Nₜ÷2:-1:2, i]
                sf[c.box.Nₜ÷2+2:end, i] .= sf[c.box.Nₜ÷2:-1:2, i] 
            end
        elseif c.seed == "cn"
            # Allocate empty arrays to hold lax pair generating functions
            rf = Array{Complex{Float64}}(undef, c.box.Nₜ, c.box.Nₓ)
            sf = Array{Complex{Float64}}(undef, c.box.Nₜ, c.box.Nₓ)
            # Set up x in the proper way
            # We currently do not suppport t shifts for this seed
            x = c.box.x' .- c.xₛ[p]

            # Set up a funnction for the ODE
            function cn_ab!(du,u,p,t)
                du[1] = +im*p[1].*u[1] .+ im.*u[2]*sqrt(p[2]).*cn.(t,real(p[2]))
                du[2] = -im*p[1].*u[2] .+ im.*u[1]*sqrt(p[2]).*cn.(t,real(p[2]))
            end
            # And the initial condition
            function ab_cn_t0(x)
                A = +c.χ[p] .+ 0.5.*c.Ω[p]*c.λ[p].*(x - c.xₛ[p]) .- π/4
                B = -c.χ[p] .+ 0.5.*c.Ω[p]*c.λ[p].*(x - c.xₛ[p]) .- π/4
                return [2*im*sin(A); 2*cos(B)]
            end

            # Prepare the ODE intnegrator
            tspan = (0.0,abs(minimum(c.box.t)))
            u0 = ab_cn_t0(c.box.x[1])
            prob = ODEProblem(cn_ab!, u0, tspan, [c.λ[p], c.m])
            dt = c.box.t[2] - c.box.t[1]
            integrator = init(prob, Tsit5(), saveat=abs.(c.box.t[c.box.Nₜ÷2+1:-1:1]), dt=dt, adaptive=false)

            # Loop over x and solve the ODE for each initial condition
            # Should do this using the ensemble interface in the future
            for i in eachindex(c.box.x)
                # Set up the new initial condition
                u0 = ab_cn_t0(c.box.x[i])
                reinit!(integrator, u0)
                # Solve
                DiffEqBase.solve!(integrator)
                # Save results
                rf[c.box.Nₜ÷2+1:-1:1, i] .= integrator.sol[1, :].*exp.(+im*(c.box.x[i])./4*(2*c.m-1))
                sf[c.box.Nₜ÷2+1:-1:1, i] .= integrator.sol[2, :].*exp.(-im*(c.box.x[i])./4*(2*c.m-1))
                # Mirror x -> x
                rf[c.box.Nₜ÷2+2:end, i] .= rf[c.box.Nₜ÷2:-1:2, i]
                sf[c.box.Nₜ÷2+2:end, i] .= sf[c.box.Nₜ÷2:-1:2, i] 
            end
        end

        # Compute ψ₁ when p = 1
        if  p == 1
            ψₜ[:,:,n+1] .= ψₜ[:,:,1] .+ (2*(c.λ[n]' - c.λ[n]).*sf.*conj.(rf)) ./
                              (abs2.(rf) + abs2.(sf))
        end
    else
        # Recursion when n != 1
        # Calculate r_{n_1, 1} and r_{n-1, p} and similarly for s
        r1, s1 = calc_rs(c, n-1, 1, ψₜ)
        r2, s2 = calc_rs(c, n-1, p+1, ψₜ)

        # Apply the DT equations
        rf = ((c.λ[n-1]'  - c.λ[n-1] ) .* conj.(s1) .* r1 .* s2 +
              (c.λ[p+n-1] - c.λ[n-1] ) .* abs2.(r1) .* r2       +
              (c.λ[p+n-1] - c.λ[n-1]') .* abs2.(s1) .* r2        ) ./
              (abs.(r1).^2 + abs.(s1).^2);
        sf = ((c.λ[n-1]'  - c.λ[n-1] ) .*       s1   .*conj(r1) .* r2 +
              (c.λ[p+n-1] - c.λ[n-1] ) .*  abs2.(s1) .*     s2        +
              (c.λ[p+n-1] - c.λ[n-1]') .*  abs2.(r1) .*     s2         ) ./
              (abs.(r1).^2 + abs.(s1).^2);
        if p == 1
            # Compute ψₙ when p = 1
            ψₜ[:,:,n+1] .= ψₜ[:,:,n] + (2*(c.λ[n]' - c.λ[n]).*sf.*conj.(rf)) ./
                                      (abs2.(rf) + abs2.(sf))
        end
    end

    return rf, sf
end