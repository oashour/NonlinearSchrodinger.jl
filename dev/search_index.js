var documenterSearchIndex = {"docs":
[{"location":"library/#Library-Documentation","page":"Library Documentation","title":"Library Documentation","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"Documentation for NLSS.jl's library.","category":"page"},{"location":"library/#Contents","page":"Library Documentation","title":"Contents","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"Pages = [\"library.md\"]","category":"page"},{"location":"library/#Index","page":"Library Documentation","title":"Index","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"Pages = [\"library.md\"]","category":"page"},{"location":"library/#Simulation","page":"Library Documentation","title":"Simulation","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"Modules = [NLSS.Simulation]\nOrder   = [:type, :function]","category":"page"},{"location":"library/#NLSS.Simulation.T₁ʰ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₁ʰ","text":"T₁ʰ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic second order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₁ˢ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₁ˢ","text":"T₁ˢ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic second order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₂ʰ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₂ʰ","text":"T₂ʰ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic second order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₂ˢ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₂ˢ","text":"T₂ˢ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic second order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₄ʰ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₄ʰ","text":"T₄ʰ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic fourth order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!, T2\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₄ˢ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₄ˢ","text":"T₄ˢ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic fourth order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!, T2\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₆ˢ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₆ˢ","text":"T₆ˢ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic sixth order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!, T₄ˢ\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.T₈ˢ-Tuple{Any,Any,Any}","page":"Library Documentation","title":"NLSS.Simulation.T₈ˢ","text":"T₈ˢ(ψ, ω, dx, F)\n\nCompute ψ', i.e. ψ advanced a step dx forward using a symplectic eighth order integrator. ψ' is defined on an FFT grid with frequencies ω using an FFT plan F. Do not call this explicitly and use solve! instead.\n\nSee also: solve!, T₆ˢ\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.print-Tuple{NLSS.Simulation.Sim}","page":"Library Documentation","title":"NLSS.Simulation.print","text":"function print(sim::Sim)\n\nPrints information about the Sim instance sim\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.solve!-Tuple{NLSS.Simulation.Sim}","page":"Library Documentation","title":"NLSS.Simulation.solve!","text":"solve!(sim::Simulation)\n\nSolves the Simulation object sim using the techniques its attributes specify.\n\nSee also: init_sim, NLSS.Plotter.plot_ψ\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Simulation.ψ₀_periodic-Tuple{Array,NLSS.Utilities.Box,Any}","page":"Library Documentation","title":"NLSS.Simulation.ψ₀_periodic","text":"function ψ₀_periodic(coeff, box::SimBox, params::SimParamseters; phase=0)\n\nComputes an initial wavefunction for the SimBox box, with fundamental frequency sim.Ω and coefficients A_1n = coeff and an overall phase exp(i phase t), i.e. of the form:\n\nTODO: insert latex form here\n\nSee also: init_sim\n\n\n\n\n\n","category":"method"},{"location":"library/#Calculation","page":"Library Documentation","title":"Calculation","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"","category":"page"},{"location":"library/#Plotter","page":"Library Documentation","title":"Plotter","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"Modules = [NLSS.Plotter]\nOrder   = [:function, :type]","category":"page"},{"location":"library/#NLSS.Plotter.plot_IoM","page":"Library Documentation","title":"NLSS.Plotter.plot_IoM","text":"function plot_IoM(sim::Sim, x_res::Int64 = 500)\n\nPlots the integrals of motion for a Sim object sim with a resolution x_res points in the x-direction. Produces plots of the energy, kinetic energy, potential energy, energy error, particle number and momentum\n\nSee also: Simulation.compute_CoM!\n\n\n\n\n\n","category":"function"},{"location":"library/#NLSS.Plotter.plot_ψ-Tuple{Any}","page":"Library Documentation","title":"NLSS.Plotter.plot_ψ","text":"plot_ψ(sim::Sim; mode = \"density\", power=1, x_res=500, t_res=512)\n\nPlot sim.abs.(ψ).^power on a grid of size approximately t_resxx_res The function can plot either a heatmap (mode = \"density\") or a 3D surface (mode = \"surface\")\n\nTODO: Insert pictures here\n\nSee also: Simulation.solve!\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Plotter.plot_ψ̃-Tuple{Any}","page":"Library Documentation","title":"NLSS.Plotter.plot_ψ̃","text":"plot_ψ̃(sim::Sim; mode = \"density\", power=1, x_res=500, t_res=512)\n\nPlot log(sim.abs.(ψ)) on a grid of size approximately ω_resxx_res if mode = density, or with n_lines total, skipping skip between each  line if mode = lines. The latter uses x_res points as well.\n\nTODO: Insert pictures here\n\nSee also: Simulation.compute_spectrum!\n\n\n\n\n\n","category":"method"},{"location":"library/#Utilities","page":"Library Documentation","title":"Utilities","text":"","category":"section"},{"location":"library/","page":"Library Documentation","title":"Library Documentation","text":"Modules = [NLSS.Utilities]\nOrder   = [:type, :function]","category":"page"},{"location":"library/#NLSS.Utilities.compute_IoM!-Tuple{Any}","page":"Library Documentation","title":"NLSS.Utilities.compute_IoM!","text":"function compute_IoM!(obj)\n\nComputes the integrals of motion of obj.ψ and saves them in respective fields of obj. obj can be a Sim or Calc object\n\nSee also: NLSS.Plotter.plot_CoM\n\n\n\n\n\n","category":"method"},{"location":"library/#NLSS.Utilities.compute_spectrum!-Tuple{Any}","page":"Library Documentation","title":"NLSS.Utilities.compute_spectrum!","text":"function compute_spectrum!(obj)\n\nComputes the normalized spectrum of obj.ψ with the center frequency shifted to the center and saves it in obj.ψ̃. obj can be a Sim or Calc object\n\nSee also: NLSS.Plotter.plot_ψ̃\n\n\n\n\n\n","category":"method"},{"location":"man/theory/#Theory","page":"Theory","title":"Theory","text":"","category":"section"},{"location":"man/theory/#Cubic-NLS","page":"Theory","title":"Cubic NLS","text":"","category":"section"},{"location":"man/theory/","page":"Theory","title":"Theory","text":"The cubic nonlinear Schrodinger equation is given by:","category":"page"},{"location":"man/theory/","page":"Theory","title":"Theory","text":"i fracpartial psipartial x + frac12 fracpartial ^2 psipartial t^2 + psi psi^2 = 0 ","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = NLSS","category":"page"},{"location":"#NLSS","page":"Home","title":"NLSS","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A suite of tools generator for the Nonlinear Schrodinger hierarchy.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for studying the nonlinear Schrodinger hierarchy's numerical and analytical solutions.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The current features are currently available or are a work in progress:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Solve the cubic nonlinear Schrodinger equation numerically using a variety of symplectic and RKN algorithms.\nCompute the integrals of motion (energy, momentum, and particle number)\nCompute the Darboux Transformation to study complicated analytical solutions of the full hierarchy\nProduce publication quality graphics using the GR backend of Plots.jl\nExport to HDF5","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Theory page provides an introduction to the theoretical backpinnings of this package to help you get started.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Some examples can be found on the Examples page.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of documented functions and types.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"man/theory.md\",\n]\nDepth = 1","category":"page"},{"location":"#Library-Outline","page":"Home","title":"Library Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"library.md\"]","category":"page"},{"location":"#main-index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"library.md\"]","category":"page"}]
}
