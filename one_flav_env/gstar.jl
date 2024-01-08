module gstar

export gstar_interp

using Interpolations
import DelimitedFiles: readdlm

# load data that was plot-grabbed from a simple plot of g* vs T
data = readdlm("$(@__DIR__)/gstar_vs_T[MeV].csv", ',', Float64)
temps = reverse(data[:,1]) .* 1e-3
gvals = reverse(data[:,2])

"""
    gstar_interp(T)

Returns the effective degrees of freedom g* for a given temperature.

T : temperature in GeV
"""
gstar_interp = linear_interpolation(temps, gvals, extrapolation_bc=Line())

end