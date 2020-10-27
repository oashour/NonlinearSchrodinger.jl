using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger(right_justify=120))
using ProgressLogging

include("test2.jl")
@progress for i=1:100
    if i == 50
        @info "Middle of computation" i
    elseif i == 70
        println("Normal output does not interfere with progress bars")
    end
    sleep(0.01)
end
@info "Done"