module LFDM_Boltz

export lifetime, lifetimes

using Pkg
Pkg.activate((@__DIR__))

using QuadGK
using DifferentialEquations
using Interpolations
using Roots
import SpecialFunctions: besselk
import LinearAlgebra: dot
import PhysicalConstants.CODATA2018: G, FineStructureConstant
using NaturallyUnitful

include((@__DIR__)*"/"*"gstar.jl")
using .gstar

include((@__DIR__)*"/"*"one_flav_model.jl")
import .OneFlavor: LFDM, Γ_φ_to_χe, σ_AA_to_φφ, σ_Ae_to_φχ, σ_ee_to_φφ, 
        σ_ZZ_to_φφ, Y_eq_ratio

# Physical Constants:
const αEM = float(FineStructureConstant) # Fine structure constant
const Grav = ustrip(uconvert(u"GeV^-2", natural(float(G)))) # Gravitational constant

const T0 = 2.348e-13 # GeV, present day temperature

const s0 = ustrip(uconvert(u"GeV^3", natural(2891.2*u"cm^-3"))) # GeV^3, 
                                                # present day entropy density

const Mstar_no_s = 3*(2.131e-42)^2/(2*T0^3*8*π*Grav) # Collection of constants
                                        # for calculating the relic abundance

const lifetime_conv = ustrip(unnatural(u"s", 1*u"GeV^-1")) # conversion from
                                            # GeV^-1 to s, used for lifetime

const h_consts = sqrt(4*π^3*Grav/45) # the constants in the hubble constant

s_factor(T) = (T0^3*2*π^2*gstar_interp(T)/45)/s0 # factor that accounts for the
                                # difference in using Y = n / T^3 and Y = n / s

Mstar_s(T) = Mstar_no_s * s_factor(T)

function arctanh(x)
    .5*(log(1+x) - log(1-x))
end

"""
    Y_eq(x, g)

Gives the equilibrium yield Y_eq = n_eq / T^3 where n_eq is the equilibrium
number density for a particle X with g spin degrees of freedom.

x : m_X / T, mass of X divided by temperature
g : # of spin degrees of freedom of X
"""
function Y_eq(x, g)
    g/(2*π^2) * x^2 * besselk(2,x)
end

function bkf(u, x)
    t1 = sqrt(x/π)
    t2 = -(57+4*u)/(16*sqrt(x*π))
    t3 = (4185+24*u*(17+2*u)) / (512 * sqrt(π*x^3))
    t4 = -5*(24135+2724*u+432*u^2+64*u^3) / (8192*sqrt(π*x^5))
    return t1+t2+t3+t4
end

function bk2_ratio(x, δ)
    if x < 1e-2
        return δ^2*(1-(1-δ^2)*x^2/4)
    elseif x < 1e2/δ
        return besselk(2,x)/besselk(2,δ*x)
    else
        return exp(-(1-δ)*x)*(sqrt(δ)-15*(1-δ)/(8*x*sqrt(δ)))
    end
end

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

# 2 to 2 cross-sections:
σ38(α, model::LFDM) = 3*4*σ_ee_to_φφ(α, model) + 4*σ_AA_to_φφ(α, model)
# the factors of 4 are for the initial spin dofs (2 for each initial particle)
# the factor of 3 is for 3 different leptons all with the same xsec

"""
    C22(σ, x, model)

The 2 to 2 collision term that enters into the Boltzmann equation.

σ : dimensionless cross-section (σ = xsec * mφ^2)
x : mφ/T
model : an instance of LFDM
"""
function C22(σ, x, model::LFDM)
    redh = reduced_H(model.mφ/x)
    if x < 1e-2
        return σ(1e10, model)*1e10/(2*model.mφ*redh*model.gφ^2)
    elseif x < 1e2
        pref = 1 / (2*x*model.mφ*redh*model.gφ^2*besselk(2,x)^2)
        I = quadgk(α -> α^(3/2)*σ(α, model)*besselk(1,x*sqrt(α)), 4, Inf)[1]
        return pref*I
    else
        pref = 1 / (x^2*model.mφ*redh*model.gφ^2)
        I = quadgk(u -> exp(-u)*(u/x+2)^4*σ((u/x+2)^2, model)*bkf(u, x), 0, Inf)[1]
        return pref*I
    end
end

"""
    C12(Γ, x, model)

The 1 to 2 collision term that enters into the Boltzmann equation.

Γ : dimension*full* decay rate [GeV]
x : mφ/T
model : an instance of LFDM
"""
function C12(Γ, x, model::LFDM)
    redh = reduced_H(model.mφ/x)
    if x < 1e-2
        return x^2*Γ(model)/(2*model.mφ^2*redh)
    elseif x < 1e2
        return x*Γ(model)*besselk(1,x)/(model.mφ^2*redh*besselk(2,x))
    else
        return x*Γ(model)/(model.mφ^2*redh) * (1-3/(2*x))
    end
end

# the Boltzmann equation:
function diffeq(du, u, p, x)
    yeq_φ = Y_eq(x, p.gφ)
    yeq_rat = Y_eq_ratio(x, p)

    decay = 3*C12(Γ_φ_to_χe, x, p)*(exp(-u[1]) - u[2]*yeq_rat)

    du[1] = exp(u[1])*(C22(σ38, x, p)*(exp(-2*u[1])-yeq_φ^2) + decay)
    du[2] = decay
end

# the Boltzmann equation for just φ, assuming no decays to χ:
function f_φ(y, p, x)
    yeq_φ = Y_eq(x, p.gφ)
    exp(y)*C22(σ38, x, p)*(exp(-2*y)-yeq_φ^2)
end

# the Boltzmann equation for just χ, assuming φ remains in equilibrium:
function f_χ(y, p, x)
    yeq_χ = Y_eq(x*p.mχ/p.mφ, p.gχ)
    yeq_rat = Y_eq_ratio(x, p)
    return 3*yeq_rat*C12(Γ_φ_to_χe, x, p)*(yeq_χ-y)
end

"""
    find_yφf(mphi; x0=1e-4, xf=1e8)

Calculates the final φ yield, assuming φ does not decay to χ
"""
function find_yφf(mphi; x0=1e-4, xf=1e8)
    M = LFDM(0.0, mphi, 0.0)

    z0 = -log(Y_eq(x0, M.gφ))

    prob = ODEProblem(f_φ, z0, (x0, xf), M)
    sol = solve(prob, Rodas4P())

    return exp(-sol(xf))
end

# integral in yχf
function yχf_I(mφ)
    quadgk(x -> x^3*besselk(1,x)/sqrt(gstar_interp(mφ/x)), 0, Inf)[1]
end

"""
    find_yχf(λ, mφ, mχ; gφ=1)

Calculates the final χ yield, assuming φ remains in equilibrium
"""
function find_yχf(λ, mφ, mχ; gφ=1)
    mm = LFDM(λ, mφ, gφ, mχ, 2)
    I = yχf_I(mφ)
    return I * gφ*3*Γ_φ_to_χe(mm) / (2*π^2*mφ^2*h_consts)
end

"""
    find_yχf(I, λ, mφ, mχ; gφ=1)

When called with the result of the integral as the first argument, skips
re-calculating the integral (useful because the integral only depends on mφ)
"""
function find_yχf(I, λ, mφ, mχ; gφ=1)
    mm = LFDM(λ, mφ, gφ, mχ, 2)

    return I * gφ*3*Γ_φ_to_χe(mm) / (2*π^2*mφ^2*h_consts)
end

"""
    find_λ(mφ, mχ)

Calculates the value of the Yukawa coupling λ that gives the correct final relic
abundance according to the approximation Y_χ^∞ = Y_χ^FI + Y_φ^FO for a single
mφ and mχ
"""
function find_λ(mφ, mχ)
    yφfo = find_yφf(mφ)
    yχfi = find_yχf(1.0, mφ, mχ)

    return sqrt( (Mstar_s(mφ/25)/mχ - yφfo) / yχfi )
end

"""
    find_λs(Mφ, Mχ)

Calculates the value of the Yukawa coupling λ that gives the correct final relic
abundance according to the approximation Y_χ^∞ = Y_χ^FI + Y_φ^FO for a vector
of mφs and mχs, and returns a matrix.
"""
function find_λs(Mφ, Mχ)
    result = zeros(typeof(Mφ[1]), (length(Mφ), length(Mχ)))
    Threads.@threads for i in eachindex(Mφ)
        yφfo = find_yφf(Mφ[i])
        yχI = yχf_I(Mφ[i])

        for j in eachindex(Mχ)
            if (Mχ[j] < Mstar_s(Mφ[i]/25)/yφfo) & (Mχ[j] < Mφ[i])
                yχfi = find_yχf(yχI, 1.0, Mφ[i], Mχ[j])
                result[i,j] = sqrt( (Mstar_s(Mφ[i]/25)/Mχ[j] - yφfo) / yχfi )
            else
                result[i,j] = NaN
            end
        end
    end
    return result
end

"""
    lifetime(λ, mφ, mχ; gφ=1, gχ=2)

Calculates the lifetime of φ given values of the Yukawa coupling λ and the two
masses.
"""
function lifetime(λ, mφ, mχ; gφ=1, gχ=2)
    mm = LFDM(λ, mφ, gφ, mχ, gχ)
    decay = 3*Γ_φ_to_χe(mm)
    return lifetime_conv/decay
end

"""
    mχ_life(τφ, mφ; gphi=1)

Calculates the mχ value that corresponds to a given lifetime and mass of φ,
such that the correct relic abundance is reproduced
"""
function mχ_life(τφ, mφ; gphi=1)
    I = yχf_I(mφ)
    γχ = I*gphi/(2*π^2*mφ^2*h_consts)

    τ = τφ/lifetime_conv
    Mstar_s(mφ/25) / (γχ/τ + find_yφf(mφ; gphi=gphi))
end

"""
    lifetime(mφ, mχ)

When called without a value for the Yukawa coupling λ, calculates the one that
reproduces the correct relic abundance.
"""
lifetime(mφ, mχ) = lifetime(find_λ(mφ, mχ), mφ, mχ)

"""
    lifetimes(Mφ, Mχ)

Calculates the lifetime at each mφ and mχ that corresponds to reproducing the
correct relic abundance. Mφ and Mχ can be vectors and this will return a matrix.
"""
lifetimes(Mφ, Mχ) = lifetime.(find_λs(Mφ, Mχ), Mφ, Mχ')

end