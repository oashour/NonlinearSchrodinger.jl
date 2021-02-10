# Examples

## Example 1: Cosine Wave Initial Condition
```@setup 1
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 1
λ, T, Ω = params(λ = 0.8im)

xᵣ = 0=>100
box = Box(xᵣ, T, dx=1e-3, Nₜ = 256, n_periods = 1)

coeff = [1e-4]
ψ₀, A₀ = ψ₀_periodic(coeff, box, Ω)

sim = Sim(λ, box, ψ₀, T4A_TJ!)

solve!(sim)
compute_IoM!(sim)
surface(sim)
savefig("example1_psi.png") #hide
heatmap(sim, :ψ̃)
savefig("example1_psi_tilde.png") #hide
plot(sim, :IoM)
savefig("example1_IoM.png") #hide
```
![](example1_psi.png)
![](example1_psi_tilde.png)
![](example1_IoM.png)

## Example 2: Soliton Initial Condition
```@setup 2
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 2
λ = 0.75im

T = 20
xᵣ = 0=>100
box = Box(xᵣ, T, dx=1e-3, Nₜ = 256, n_periods = 1)

ψ₀ = Array{Complex{Float64}}(undef, box.Nₜ)
ψ₀ .= 2*imag(λ)./cosh.(2*imag(λ).*box.t)

sim = Sim(λ, box, ψ₀, T4A_TJ!)

solve!(sim)
surface(sim)
savefig("example2.png") # hide
```
![](example2.png)

## Exampole 3: Nonlinear Talbot Carpet (Pruning)
```@setup 3
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 3
λ, T, Ω = params(a = 0.36)

xᵣ = 0=>60
box = Box(xᵣ, T, dx=1e-4, Nₜ = 512, n_periods = 5)

coeff = [(2.7 + 4.6im)*1e-2]
ψ₀, A₀ = ψ₀_periodic(coeff, box, Ω)

sim = Sim(λ, box, ψ₀, T4A_TJ!, β=10.0)

solve!(sim)
surface(sim)
savefig("example3_psi.png") #hide
plot(sim, :ψ̃)
savefig("example3_psi_tilde.png") #hide
```
![](example3_psi.png)
![](example3_psi_tilde.png)

## Example 4: 7 Soliton Collision
```@setup 4
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 4
xᵣ = -10=>10
T = 20
seed = "0"
box = Box(xᵣ, T, Nₓ=1000, Nₜ = 1000)
λ = [-0.45 + 0.775im, -0.35 + 0.8im, -0.25 + 0.825im, 0.85im, 0.25 + 0.875im, 0.35 + 0.9im, 0.45 + 0.925im]
xₛ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

calc = Calc(λ, tₛ, xₛ, seed, box) 

solve!(calc)
heatmap(calc) 
savefig("example4.png") #hide
```
![](example4.png)

## Example 5: Fifth Order Maximal Intensity Breather
```@setup 5
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 5
xᵣ = -10=>10
λ₁ = 0.98im
λ, T, Ω = params(λ = λ₁)
box = Box(xᵣ, T, Nₓ=1000, Nₜ = 1000, n_periods = 3)

λ = λ_maximal(λ₁, 5) # array of 5 eigenvalues
xₛ = [0.0, 0.0, 0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0, 0.0, 0.0]

seed = "exp"
calc = Calc(λ, tₛ, xₛ, seed, box) 

solve!(calc)
surface(calc) 
savefig("example5.png") #hide
```
![](example5.png)

## Example 6: 3 Soliton Collision on a ``cn`` background
```@setup 6
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 6
xᵣ = -10=>10
T = 20

box = Box(xᵣ, T, Nₓ=1000, Nₜ = 1000, n_periods = 1)
λ = [-0.3+0.85im, 0.9im, 0.3+0.95im]
xₛ = [0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0]

seed = "cn"
calc = Calc(λ, tₛ, xₛ, seed, box, m = 0.5) 

solve!(calc)
surface(calc) 
savefig("example6.png") #hide
```
![](example6.png)

## Example 7: First Order Breather mMtched to a ``dn`` Background
```@setup 7
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 7
xᵣ = -10=>10
m = 2/5
λ = λ_given_m(m, q=4)
λ, T, Ω = params(λ = λ, m=m)
box = Box(xᵣ, T, Nₓ=1000, Nₜ = 1000, n_periods = 3)

λ = λ_maximal(λ, 1, m=m)
xₛ = [0.0]
tₛ = [0.0]

seed = "dn"
calc = Calc(λ, tₛ, xₛ, seed, box, m=m) 

solve!(calc)
heatmap(calc) 
savefig("example7.png") #hide
```
![](example7.png)

## Example 8: Darboux Transformation Initial Condition
```@setup 8
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 8
λ₁ = 0.98im
λ, T, Ω = params(λ = λ₁)

xᵣ = 0=>100
box = Box(xᵣ, T, dx=1e-3, Nₜ = 512, n_periods = 1)

λ = λ_maximal(λ₁, 5) # array of 5 eigenvalues
xₛ = [0.0, 0.0, 0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0, 0.0, 0.0]
ψ₀ = ψ₀_DT(λ, tₛ, xₛ, -10, box)

sim = Sim(λ₁, box, ψ₀, T4A_TJ!)

solve!(sim)
compute_IoM!(sim)
surface(sim) 
savefig("example8.png") #hide
```
![](example8.png)