# Darboux Transformation Examples

## Example 1: 7 Soliton Collision
```@setup 1
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 1
xᵣ = -10=>10
T = 20
seed = "0"
box = Box(xᵣ, T, Nₓ=500, Nₜ = 500)
λ = [-0.45 + 0.775im, -0.35 + 0.8im, -0.25 + 0.825im, 0.85im, 0.25 + 0.875im, 0.35 + 0.9im, 0.45 + 0.925im]
xₛ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

calc = Calc(λ, tₛ, xₛ, seed, box) 

solve!(calc)
heatmap(calc) 
savefig("example1.svg") #hide
```
![](example1.svg)

## Example 2: Fifth Order Maximal Intensity Breather
```@setup 2
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 2
xᵣ = -10=>10
λ₁ = 0.98im
λ, T, Ω = params(λ = λ₁)
box = Box(xᵣ, T, Nₓ=500, Nₜ = 500, n_periods = 3)

λ = λ_maximal(λ₁, 5) # array of 5 eigenvalues
xₛ = [0.0, 0.0, 0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0, 0.0, 0.0]

seed = "exp"
calc = Calc(λ, tₛ, xₛ, seed, box) 

solve!(calc)
surface(calc) 
savefig("example2.svg") #hide
```
![](example2.svg)

## Example 3: 3 Soliton Collision on a ``cn`` background
```@setup 3
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 3
xᵣ = -10=>10
T = 20

box = Box(xᵣ, T, Nₓ=500, Nₜ = 500, n_periods = 1)
λ = [-0.3+0.85im, 0.9im, 0.3+0.95im]
xₛ = [0.0, 0.0, 0.0]
tₛ = [0.0, 0.0, 0.0]

seed = "cn"
calc = Calc(λ, tₛ, xₛ, seed, box, m = 0.5) 

solve!(calc)
surface(calc) 
savefig("example3.svg") #hide
```
![](example3.svg)

## Example 4: First Order Breather matched to a ``dn`` Background
```@setup 4 
using NonlinearSchrodinger
using Plots
using LaTeXStrings
```

```@example 4
xᵣ = -10=>10
m = 2/5
λ = λ_given_m(m, q=4)
λ, T, Ω = params(λ = λ, m=m)
box = Box(xᵣ, T, Nₓ=500, Nₜ = 500, n_periods = 3)

λ = λ_maximal(λ, 1, m=m)
xₛ = [0.0]
tₛ = [0.0]

seed = "dn"
calc = Calc(λ, tₛ, xₛ, seed, box, m=m) 

solve!(calc)
heatmap(calc) 
savefig("example4.svg") #hide
```
![](example4.svg)