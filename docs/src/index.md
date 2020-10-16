```@meta
CurrentModule = NLSS
```
# NLSS

*A suite of tools generator for the Nonlinear Schrodinger hierarchy.*

A package for studying the nonlinear Schrodinger hierarchy's numerical and analytical solutions.

## Package Features
The current features are currently available or are a work in progress:

- Solve the cubic nonlinear Schrodinger equation numerically using a variety of symplectic and RKN algorithms.
- Compute the integrals of motion (energy, momentum, and particle number)
- Compute the Darboux Transformation to study complicated analytical solutions of the full hierarchy
- Produce publication quality graphics using the GR backend of Plots.jl
- Export to HDF5

The [Theory](@ref) page provides an introduction to the theoretical backpinnings of this package to help you get started.

Some examples can be found on the [Examples](@ref) page.

See the [Index](@ref main-index) for the complete list of documented functions and types.
```@contents
Pages = [
    "man/theory.md",
]
Depth = 1
```

## Library Outline

```@contents
Pages = ["lib/public.md", "lib/private.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```
