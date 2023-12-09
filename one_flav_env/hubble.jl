module Hubble

export hubble, reduced_H

import PhysicalConstants.CODATA2018: G
using NaturallyUnitful

include((@__DIR__)*"/"*"gstar.jl")
using .gstar

const Grav = ustrip(uconvert(u"GeV^-2", natural(float(G)))) # Gravitational constant

"""
    hubble(T)

Gives the Hubble constant H in GeV.

T : Temperature [GeV]
"""
function hubble(T)
    T^2 * sqrt(4*π^3*Grav*gstar_interp(T)/45)
end

"""
    reduced_H(T)

Gives H / T^2 in GeV^-1. For use in Boltzmann equations.

T : Temperature [GeV]
"""
function reduced_H(T)
    sqrt(4*π^3*Grav*gstar_interp(T)/45)
end

end