```@meta
CurrentModule = NonlinearSchrodinger
```
# NonlinearSchrodinger

*A suite of tools for Nonlinear Schrodinger Equations*

`NonlinearSchrodinger.jl` is a suite of tools for solving Nonlinear Schrodinger equations via higher-order algorithms and Darboux transformations.

## Package Features
The following features are currently available:

- Solving the cubic Nonlinear Schrodinger equation using a plethora of algortithms of order up to 8 (the number of algorithms available is always increasing!). Symplectic and Nystrom integrators are available.
- Solving the Hirota and Sasa-Satsuma equations using a combined split-step-finite-difference approach using a few different integrators. 
- Computing the integrals of motion (energy, momentum, and particle number) and their errors.
- Computing the Darboux Transformation to study complicated analytical solutions. We currently support the breather and soliton seeds for extended nonlinear Schrodinger equations of order up to 5 (including cubic NLS, Hirota, LPD, Quintic, and arbitrary combinations thereof). We also support the `cn` and `dn` seeds for the cubic NLS.
- Easy [Visualization](@ref) through `Plots.jl` recipes.
- Very simple API that allows one to compute very complicated solutions via only a few lines of code. Some examples can be found on the [Examples](@ref) page.


```@contents
Pages = [
    "man/theory.md",
]
Depth = 1
```

## Library Outline

```@contents
Pages = ["public.md"]
```

### [Index](@id main-index)

```@index
Pages = ["public.md"]
```
