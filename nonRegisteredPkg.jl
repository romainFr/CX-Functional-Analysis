Pkg.build("PlotlyJS")
Pkg.clone("https://github.com/romainFr/FluorescentSeries.jl.git")
using Blink
Blink.AtomShell.install()
