using Plots
using FFTW
using NonlinearSchrodinger

function T1A_K(ψᵢ, b)
    γ = 1
    p = 1
    q = 1

    n = 1
    μ = 0 

    #ψ = ifft(im*b.ω.*fftshift(fft(ψᵢ)))
    #ψₒ = exp.(-γ*ψᵢ.^(n-1).*b.dx.*ψ).*ψᵢ
    #ψₒ = fftshift(fft(ψₒ))
    #ψₒ = ifft(exp.(-b.dx*(im*b.ω*μ .+ im*b.ω.^3*p .+ im*b.ω.^5*q)).*ψₒ)

    ψₒ = fftshift(fft(ψᵢ))
    ψₒ = ifft(exp.(-b.dx*(im*b.ω*μ .+ im*b.ω.^3*p .+ im*b.ω.^5*q)).*ψₒ)
    ψₒ = exp.(-γ*ψᵢ.^(n-1).*b.dx.*ifft(+im*b.ω.*fftshift(fft(ψₒ)))).*ψₒ
end

xᵣ = 0=>50
box = Box(xᵣ, 200, dx=5e-2, Nₜ = 512, n_periods =3)

ψ = zeros(box.Nₓ, box.Nₜ) .+ im

ψ[1, :] = 105/169*sech.(1/2/sqrt(13)*box.t).^4
for i = 1:box.Nₓ-1
   ψ[i+1, :] .=  T1A_K(ψ[i,:], box)
end
heatmap(abs.(ψ))