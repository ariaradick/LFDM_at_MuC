module OneFlavor

export LFDM, Γ_φ_to_χe, σ_AA_to_φφ, σ_Ae_to_φχ, σ_ee_to_φφ, σ_ZZ_to_φφ,
        Y_eq, Y_eq_ratio, C12, C22

using QuadGK
import SpecialFunctions: besselk
import PhysicalConstants.CODATA2018: FineStructureConstant
using NaturallyUnitful

include((@__DIR__)*"/"*"hubble.jl")
using .Hubble

const αEM = float(FineStructureConstant)
const sW2 = .23121 # sine of weak-mixing angle squared
const cW2 = 1 - sW2
const mZ = 91.1876 # GeV
const mW = 80.370 # GeV
# const gc20 = sqrt(4*π*αEM) * (sW2-.5) / (sqrt(sW2*cW2))
const gc20 = sqrt(4*π*αEM*sW2/cW2)

function arctanh(x)
    .5*(log(1+x) - log(1-x))
end

function arccoth(x)
    .5*(log(x+1) - log(x-1))
end

"""
Struct that includes relevant parameters for the one flavor model.
Most observables take this as an input.

λ : Yukawa coupling of operator φ \bar{l} χ

mφ : Mass of new charged scalar φ

gφ : spin degrees of freedom for φ, should be 1

mχ : Mass of new fermion χ

gχ : spin degrees of freedom for χ, should be 2

h : Hubble evaluated at mφ
"""
struct LFDM{T <: Real}
    λ::T
    mφ::T
    gφ::T
    mχ::T
    gχ::T
    h::T
end

function LFDM(λ::Real, mφ::Real, mχ::Real)
    hub = hubble(mφ)
    LFDM(promote(λ, mφ, 1, mχ, 2, hub)...)
end

function LFDM(λ::Real, mφ::Real, gφ::Real, mχ::Real, gχ::Real)
    hub = hubble(mφ)
    LFDM(promote(λ, mφ, gφ, mχ, gχ, hub)...)
end

"""
    Γ_φ_to_χe(model::LFDM)

Decay rate of φ to χ + l (no distinction between leptons in this model). 

Units will be same as model.mφ
"""
function Γ_φ_to_χe(model::LFDM)
    model.λ^2*model.mφ*(1-model.mχ^2/model.mφ^2)^2 / (16*π)
end

# note: cross-sections are σ = xsec * mφ^2

"""
    σ_ee_to_φφ(α, model::LFDM)

Dimensionless cross-section for l+ l- to φ+ φ- (σ = xsec * mφ^2)
Input: α = s / mφ^2
"""
function σ_ee_to_φφ(α, model::LFDM)
    η = sqrt(1-4/α)
    z = mZ / model.mφ
    inside = 5*α^2-12*α*cW2*(α-2*sW2*(α-z^2)) + 8*cW2^2*(α-2*sW2*(α-z^2))^2
    pref = αEM^2*π/(96*sW2^2*cW2^2)
    return pref*inside*η^3 / (α*(α-z^2)^2)
end

"""
    σ_AA_to_φφ(α, model::LFDM)

Dimensionless cross-section for diphoton to φ+ φ- (σ = xsec * mφ^2)
Input: α = s / mφ^2
Takes a model as input for simplification of Boltzmann code, although not used.
"""
function σ_AA_to_φφ(α, model::LFDM)
    β = 4/α
    η = sqrt(1-4/α)
    return 2*π*αEM^2/α * (η*(2+β)-β*(1+η^2)*arctanh(η))
end

"""
    σ_ZZ_to_φφ(α, model::LFDM)

Dimensionless cross-section for di-Z to φ+ φ- (σ = xsec * mφ^2)
Input: 
α = s / mφ^2 : The dimensionless center-of-mass energy squared
model : An instance of LFDM
"""
function σ_ZZ_to_φφ(α, model::LFDM; g_Zφφ=gc20)
    z = mZ / model.mφ
    η = sqrt(1-4/α)

    sw = sqrt(sW2)
    cw = sqrt(cW2)

    pref = -1 / (144*π*α*z^4*(α-4*z^2))

    t1 = -16*π^2*αEM^2*η*sW2^2*sqrt(α*(α-4*z^2))*(α^2+12*z^4-4*α*z^2) / cW2^2
    t2 = (8*π*αEM*g_Zφφ^2*sw^2*(η*sqrt(α*(α-4*z^2))*(α-2*z^2)^2+4*z^4*
        (α+2*z^2-8)*arccoth((α-2*z^2)/sqrt(α*η^2*(α-4*z^2)))))/cw^2
    t3 = (g_Zφφ^4 / ((2*z^2-α)*(α+z^4-4*z^2)))*(12*z^12*
        log((α-η*sqrt(α*(α-4*z^2))-2*z^2)/(α+η*sqrt(α*(α-4*z^2))-2*z^2))+
        8*z^4*(α*((α-16)*α+16)+4*(α+5)*z^6+((α-35)*α-16)*z^4+8*(9*α-8)*z^2)*
        arccoth((α-2*z^2)/sqrt(α*η^2*(α-4*z^2)))-η*(2*z^2-α)*sqrt(α*(α-4*z^2))*
        (α^3+6*z^8-4*(α+8)*z^6+(α*(α+20)+32)*z^4-8*α^2*z^2))
    
    return pref*(t1+t2+t3)
end

"""
    σ_WW_to_φφ(α, model::LFDM)

Dimensionless cross-section for W+ W- to φ+ φ- (σ = xsec * mφ^2)
Input: 
α = s / mφ^2 : The dimensionless center-of-mass energy squared
model : An instance of LFDM
"""
function σ_WW_to_φφ(α, model::LFDM; g_Zφφ=gc20)
    z = mZ / model.mφ
    w = mW / model.mφ
    η = sqrt(1 - 4/α)

    pref = αEM*η^3*sqrt(α*(α-4*w^2)) / (432*α^2*sW2*w^4*(α-z^2)^2)
    t1 = α^2 + 12*w^4 + 20*α*w^2
    t2 = (α*sqrt(cW2)*g_Zφφ - 2*sqrt(π*αEM*sW2)*(α-z^2))^2

    return pref*t1*t2
end

"""
    σ_Ae_to_φχ(α, model::LFDM)

Dimensionless cross-section for photon lepton(+/-) to χ φ(+/-) (σ = xsec * mφ^2)
Input: α = s / mφ^2
"""
function σ_Ae_to_φχ(α, model::LFDM)
    δ = model.mχ/model.mφ
    sqrt_term = sqrt((α-1)^2-2*(1+α)*δ^2+δ^4)
    pref = αEM*model.λ^2 / (16*α^3)
    t1 = sqrt_term*(5*α+δ^2-1)
    t2 = 8*(2*α-(1-δ^2)^2)*arccoth( (1+α-δ^2) / sqrt_term )
    return pref*(t1-t2)
end

# Model-dependent Boltzmann equation stuff:
function bk2_ratio(x, δ)
    if x < 1e-2
        return δ^2*(1-(1-δ^2)*x^2/4)
    elseif x < 1e2/δ
        return besselk(2,x)/besselk(2,δ*x)
    else
        return exp(-(1-δ)*x)*(sqrt(δ)-15*(1-δ)/(8*x*sqrt(δ)))
    end
end

function Y_eq_ratio(x, model::LFDM)
    (model.gφ/model.gχ)*(model.mφ/model.mχ)^2*bk2_ratio(x, model.mχ/model.mφ)
end

end