module NonlinearSchrodinger

using FFTW
using OrdinaryDiffEq, DiffEqOperators
using Memoization
using LaTeXStrings  
using Elliptic
using Elliptic.Jacobi
using RecipesBase
using PolynomialRoots

export compute_IoM!, compute_spectrum!
export print
export solve!
export ψ₀_periodic, ψ₀_DT
export λ_maximal, λ_given_m, λ_given_f
export params
export PHF, find_peaks
#export solve_dt!

export Density_ψ

export Sim, Box, Calc
export T1A!, T1B!, T2A!, T2B!
export T4A_TJ!, T4B_TJ!, T6A_TJ!, T6B_TJ!, T8A_TJ!, T8B_TJ!
export T4A_SF!, T4B_SF!, T6A_SF!, T6B_SF!, T8A_SF!, T8B_SF!
export T4A_CMP!, T4B_CMP!, T6A_CMP!, T6B_CMP!, T8A_CMP!, T8B_CMP!
export T6A_Ss14!, T6B_Ss14!, T6A_Ys7!, T6B_Ys7!, T8A_Ss15!, T8B_Ss15!
export T1A_H!, T2A_H!
export T1A_SS!

include("Types.jl")
include("Simulation.jl")
include("Utilities.jl")
include("Recipes.jl")
include("Darboux.jl")


end #module
