module NLSS

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger(right_justify=120))

include("Utilities.jl")
include("Simulation.jl")
include("Plotter.jl")

end #module
