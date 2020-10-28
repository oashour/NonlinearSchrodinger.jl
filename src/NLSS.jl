module NLSS

using FFTW
using DifferentialEquations, DiffEqOperators, NumericalIntegration
using Memoization
using Plots; gr()
using Plots.PlotMeasures, LaTeXStrings  
using JLD

# Logging Stuff
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using ProgressLogging

export compute_IoM!, compute_spectrum!
export plot_ψ, plot_ψ̃, plot_IoM
export params, print, save, load
export solve!
export ψ₀_periodic
export AlgorithmType, AlgorithmVariant

export Sim, Box, Algorithm
export T1A!, T1B!, T2A!, T2B!
export T4A_TJ!, T4B_TJ!, T6A_TJ!, T6B_TJ!, T8A_TJ!, T8B_TJ!
export T4A_SF!, T4B_SF!, T6A_SF!, T6B_SF!, T8A_SF!, T8B_SF!
export T4A_N!, T4B_N!, T6A_N!, T6B_N!, T8A_N!, T8B_N!
export T6A_OP!, T6B_OP!, T8A_OP!, T8B_OP!
export T1A_H, T2A_H, T4A_TJ_H

include("Types.jl")
include("Simulation.jl")
include("Utilities.jl")
include("Plotter.jl")

end #module
