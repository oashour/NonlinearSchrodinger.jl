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

export Sim, Box

include("Types.jl")
include("Simulation.jl")
include("Utilities.jl")
include("Plotter.jl")

end #module
