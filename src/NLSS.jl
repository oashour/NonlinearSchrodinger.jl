module NLSS

using FFTW
using Plots; gr()
Plots.GRBackend()
using LaTeXStrings
using ProgressMeter

include("sim_structs.jl")
include("sim_utils.jl")
include("sim_solvers.jl")
include("plot.jl")

end #module
