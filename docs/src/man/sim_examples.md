# Simulation Examples

## Example 1: Cosine Wave Initial Condition
```@setup 1
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 1
λ, T, Ω = params(λ = 0.8im)

xᵣ = 0=>20
box = Box(xᵣ, T, dx=5e-3, Nₜ = 256, n_periods = 1)

coeff = [1e-4]
ψ₀, A₀ = ψ₀_periodic(coeff, box, Ω)

sim = Sim(λ, box, ψ₀, T4A_TJ!)

solve!(sim)
compute_IoM!(sim)
surface(sim)
savefig("example1_psi.svg") #hide
heatmap(sim, :ψ̃)
savefig("example1_psi_tilde.svg") #hide
plot(sim, :IoM)
savefig("example1_IoM.svg") #hide
```
![](example1_psi.svg)
![](example1_psi_tilde.svg)
![](example1_IoM.svg)

## Example 2: Soliton Initial Condition
```@setup 2
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 2
λ = 0.75im

T = 20
xᵣ = 0=>20
box = Box(xᵣ, T, dx=1e-2, Nₜ = 256, n_periods = 1)

ψ₀ = Array{Complex{Float64}}(undef, box.Nₜ)
ψ₀ .= 2*imag(λ)./cosh.(2*imag(λ).*box.t)

sim = Sim(λ, box, ψ₀, T4A_TJ!)

solve!(sim)
surface(sim)
savefig("example2.svg") # hide
```
![](example2.svg)

## Example 3: Nonlinear Talbot Carpet (Pruning)
```@setup 3
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 3
λ, T, Ω = params(a = 0.36)

xᵣ = 0=>60
box = Box(xᵣ, T, dx=1e-2, Nₜ = 512, n_periods = 5)

coeff = [(2.7 + 4.6im)*1e-2]
ψ₀, A₀ = ψ₀_periodic(coeff, box, Ω)

sim = Sim(λ, box, ψ₀, T4A_TJ!, β=10.0)

solve!(sim)
heatmap(sim)
savefig("example3_psi.svg") #hide
plot(sim, :ψ̃)
savefig("example3_psi_tilde.svg") #hide
```
![](example3_psi.svg)
![](example3_psi_tilde.svg)

## Example 4: Darboux Transformation Initial Condition
```@setup 4
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 4
λ₁ = 0.98im
λ, T, Ω = params(λ = λ₁)

xᵣ = 0=>20
box = Box(xᵣ, T, dx=5e-3, Nₜ = 512, n_periods = 1)

λ = λ_maximal(λ₁, 5) # array of 5 eigenvalues
xₛ = [0.0, 0.0, 0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0, 0.0, 0.0]
ψ₀ = ψ₀_DT(λ, tₛ, xₛ, -10, box)

sim = Sim(λ₁, box, ψ₀, T4A_TJ!)

solve!(sim)
compute_IoM!(sim)
heatmap(sim) 
savefig("example4.svg") #hide
```
![](example4.svg)