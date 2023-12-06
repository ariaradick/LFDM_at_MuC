"""
Implementation of EWA PDFs following arXiv:[2111.02442]
"""
module MuPDF
export σ_hadronize

using HCubature

# particles:
# {0,1,2,3,4,5,6} = {g, d, u, s, c, b, t}
# {11,12,13,14,15,16} = {e-, νe, μ-, νμ, τ-, ντ}
# {21,22,23,24,25} = {g, γ, Z, W+, h}

# [g_R, g_L, g_V, g_A]
const g_γ = [1.0, 1.0, 1.0, 0.0]
const g_W = [0.0, 1.0, 0.5, -0.5]

const sw2 = .23121
const mμ = 10^3*105.6583745e-3 # GeV
const mW = 80.377 # GeV
const mZ = 91.1876 # GeV

Qf(pdgid::Int) = -(pdgid % 2) # only valid for pdgid in 11:16
T3fL(pdgid::Int) = -((pdgid % 2) - sign(pdgid)*1/2) # only valid for pdgid in 11:16

gRf_Z(pdgid) = -Qf(pdgid)*sw2
gLf_Z(pdgid) = T3fL(pdgid) - Qf(pdgid)*sw2
gVf_Z(pdgid) = .5*T3fL(pdgid) - Qf(pdgid)*sw2
gAf_Z(pdgid) = -.5*T3fL(pdgid)

# [g_R, g_L, g_V, g_A]
g_Z(pdgid) = [gRf_Z(pdgid), gLf_Z(pdgid), gVf_Z(pdgid), gAf_Z(pdgid)]

f_Vp_fL(ξ, μf2, mV, gV, gL) = gV^2*gL^2*(1-ξ)^2*log(μf2/mV^2) / (8*π^2*ξ)
f_Vm_fL(ξ, μf2, mV, gV, gL) = gV^2*gL^2*log(μf2/mV^2) / (8*π^2*ξ)
f_V0_fL(ξ, μf2, gV, gL) = gV^2*gL^2*(1-ξ) / (4*π^2*ξ)

f_Vp(ξ, μf2, mV, g) = .5*(f_Vp_fL(ξ, μf2, mV, g[3], g[2]) + f_Vm_fL(ξ, μf2, mV, g[3], g[1]))
f_Vm(ξ, μf2, mV, g) = .5*(f_Vm_fL(ξ, μf2, mV, g[3], g[2]) + f_Vp_fL(ξ, μf2, mV, g[3], g[1]))
f_V0(ξ, μf2, g) = .5*(f_V0_fL(ξ, μf2, g[3], g[2]) + f_V0_fL(ξ, μf2, g[3], g[1]))

function pdf_γ(ξ, μf2)
    return (f_Vp(ξ, μf2, mμ, g_γ) + f_Vm(ξ, μf2, mμ, g_γ))
end

function pdf_Z(ξ, μf2; pdgid=13)
    gz = g_Z(pdgid)
    return (f_Vp(ξ, μf2, mZ, gz) + f_Vm(ξ, μf2, mZ, gz) + f_V0(ξ, μf2, gz))
end

function pdf_WT(ξ, μf2)
    return (f_Vp(ξ, μf2, mW, g_W) + f_Vm(ξ, μf2, mW, g_W))
end

function pdf_WL(ξ, μf2)
    return (f_V0(ξ, μf2, g_W))
end

pdf_W(ξ, μf2) = pdf_WT(ξ, μf2) + pdf_WL(ξ, μf2)

PDFs = Dict(22 => pdf_γ, 23 => pdf_Z, 24 => pdf_W)

function lumi_varch(x, partonid, Q2pdf)
    pdf(x, Q2pdf) = x*PDFs[partonid](x, Q2pdf)
    u = x[1]
    v = x[2]
    return -2*u*pdf(exp(u*v), Q2pdf)*pdf(exp(u-u*v), Q2pdf)
end

function σ_hadronize(σfunc, α, τmin, partonid, Q2pdf; rtol=1.49e-6)
    integrand(x) = lumi_varch(x, partonid, Q2pdf) * σfunc(α*exp(x[1]))
    return hcubature(integrand, [log(τmin), 0], [0, 1], rtol=rtol)[1]
end

end