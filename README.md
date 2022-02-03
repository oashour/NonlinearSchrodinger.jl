![Logo](logo.png?raw=true "NonlinearSchrodinger.jl Logo")
# NonlinearSchrodinger.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://oashour.github.io/NonlinearSchrodinger.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://oashour.github.io/NonlinearSchrodinger.jl/dev)
[![Build Status](https://travis-ci.com/oashour/NonlinearSchrodinger.jl.svg?branch=master)](https://travis-ci.com/oashour/NonlinearSchrodinger.jl)
[![CI](https://github.com/oashour/NonlinearSchrodinger.jl/workflows/CI/badge.svg)](https://github.com/oashour/NonlinearSchrodinger.jl/actions)
[![Coverage](https://codecov.io/gh/oashour/NonlinearSchrodinger.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/oashour/NonlinearSchrodinger.jl)
[![arXiv](https://img.shields.io/badge/arXiv-2103.14469-b31b1b.svg?style=plastic?logo=arxiv)](https://arxiv.org/abs/2103.14469)

NonlinearSchrodinger.jl is a suite of tools for studying Nonlinear Schrodinger equations. The purpose of this package is to encourage the use of open source software in studying these equations as most works in this field rely on closed-source codes. This allows for reproducability, lowers the barrier for new researchers, and alleviates the need to reinvent the wheel.


## Features

1. Solving the cubic Nonlinear Schrodinger equation using a plethora of algortithms of order up to 8 (the number of algorithms available is always increasing!). Symplectic and RKN integrators are available.

2. Solving the Hirota and Sasa-Satsuma equations using a combined split-step-finite-difference approach using a few different integrators. This is a preliminary feature and is not yet fully supported. 

3. Computing the integrals of motion (energy, momentum, and particle number) and their errors.

4. Computing the Darboux Transformation to study complicated analytical solutions. We currently support the breather and soliton seeds for extended nonlinear Schrodinger equations of order up to 5 (including cubic NLS, Hirota, LPD, Quintic, and arbitrary combinations thereof). We also support the cn and dn seeds for the cubic NLS.

5. Easy visualization through Plots.jl recipes.

6. A very simple interface that allows one to compute very complicated solutions via only a few lines of code.

7. Many utilities for studying maximal intensity breather families on uniform and dnoidal backgrounds, pruning for Nonlinear talbot carpets, and breather to soliton conversion in extended NLSEs.

## Example: 7 Soliton Collision
```
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
```
![Logo](example.png?raw=true "Example: 7 Soliton Collision")

## Supporting and Citing

The paper can be found [here](https://arxiv.org/abs/2103.14469). It is currently under review.

## Logo

The logo is formed by 3 Akhmediev breathers with a = 3/8, an homage to my first paper in the field.
