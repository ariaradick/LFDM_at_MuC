# This file will calculate the dilepton (phi) kinematic variables for each point
# for the prompt (long-lived) cases and save those as "dilepton_summary.csv"
# ("diphi_summary.csv") in the respective run folder. Also organizes the rts,
# mphi, and cross-section values for each run into "run_info.csv" at the base
# of each madgraph output folder.

using Pkg
Pkg.activate((@__DIR__)*"/../one_flav_env")

using LHEF
using LorentzVectorHEP
using DataFrames
using CSV
using StaticArrays
using LinearAlgebra
using DelimitedFiles: readdlm
using Optim

MG_RESULTS_DIR = (@__DIR__)*"/../madgraph/data/"

# make sure these match what the script has output
DIR_llνν = "mumu_to_llvv/"
DIR_φφ = "mumu_to_phiphi/"
DIR_llχχ = "mumu_to_llchichi/"
DIR_VVφφ = "VV_to_phiphi/"
DIR_VVχχ = "VV_to_llchichi/"

FOLDERS = MG_RESULTS_DIR .* ["prompt/", "prompt_scan/", "LLP/", "LLP_scan/"]

XSEC_RESULTS_FILENAME = "xsec.txt"
LHE_FILENAME = "unweighted_events.lhe.gz"

function dir_runx(x)
    if x < 10
        return "Events/run_0$x/"
    else
        return "Events/run_$x/"
    end
end

lhe_v4(p) = LorentzVector(p.e, p.px, p.py, p.pz)
lhe_v3(p) = SA[p.px, p.py, p.pz]
lhe_pt(p) = SA[p.px, p.py]

gamma(lv::LorentzVector) = energy(lv)/mass(lv)
beta(lv::LorentzVector) = sqrt(1-1/gamma(lv)^2)

function trans_mass(pT1, pT2, m1, m2)
    ET1 = sqrt((pT1 ⋅ pT1) + m1^2)
    ET2 = sqrt((pT2 ⋅ pT2) + m2^2)
    return m1^2 + m2^2 + 2*(ET1*ET2 - (pT1 ⋅ pT2))
end

function MT2(particle_1, particle_2)
    m1 = particle_1.m
    m2 = particle_2.m

    pT1 = lhe_pt(particle_1)
    pT2 = lhe_pt(particle_2)
    pt_miss = -(pT1 .+ pT2)
    init_x = MArray{Tuple{2}}(pt_miss) ./ 2

    function mt_to_min(x)
        max(trans_mass(pT1, x, m1, 0), trans_mass(pT2, pt_miss-x, m2, 0))
    end

    # res = optimize(mt_to_min, init_x, LBFGS(); autodiff=:forward)
    res = optimize(mt_to_min, init_x, NelderMead())

    return sqrt(minimum(res))
end

function find_leptons(event)
    idxs = Vector{Int64}()
    for p in event.particles[3:end]
        if p.id in [11,-11,13,-13]
            append!(idxs, p.idx)
        end
    end
    return idxs
end

function check_leptons(events)
    for (i,e) in enumerate(events)
        if length(find_leptons(e)) > 2
            throw(ArgumentError("Too many leptons in event $i"))
        end
    end
end

function dilepton_summarize(event::LHEF.Event)
    l_idx = find_leptons(event)

    p1 = lhe_v4(event.particles[l_idx[1]])
    p2 = lhe_v4(event.particles[l_idx[2]])

    p1_3 = lhe_v3(event.particles[l_idx[1]])
    p2_3 = lhe_v3(event.particles[l_idx[2]])

    inv_mass = mass(p1+p2)
    pT = pt(p1+p2)
    open_angle = (p1_3 ⋅ p2_3) / (norm(p1_3)*norm(p2_3))
    deltaR = ΔR(p1, p2)
    mt2 = MT2(event.particles[l_idx[1]], event.particles[l_idx[2]])
    
    return [inv_mass, pT, open_angle, deltaR, mt2]
end

function dilepton_summarize(lhe::Vector{LHEF.Event})
    N_events = length(lhe)

    res = zeros(Float64, (N_events, 5))

    for i in 1:N_events
        res[i,:] = dilepton_summarize(lhe[i])
    end

    res = Matrix(res')

    return res
end

function dilepton_pT(event)
    l_idx = find_leptons(event)

    p1 = lhe_v4(event.particles[l_idx[1]])
    p2 = lhe_v4(event.particles[l_idx[2]])

    return pt(p1+p2)
end

function find_phis(event)
    idxs = Vector{Int64}()
    for p in event.particles[3:end]
        if p.id in [9000006, -9000006]
            append!(idxs, p.idx)
        end
    end
    return idxs
end

function phi_summarize(event::LHEF.Event)
    phi_idx = find_phis(event)

    p1 = lhe_v4(event.particles[phi_idx[1]])
    p2 = lhe_v4(event.particles[phi_idx[2]])

    gammas = gamma.([p1,p2])
    betas = beta.([p1,p2])
    betagamma = betas .* gammas

    etas = eta.([p1,p2])
    
    return Matrix(reduce(hcat, [betas, betagamma, etas])')
end

function phi_summarize(lhe::Vector{LHEF.Event})
    N_events = length(lhe)

    res = zeros(Float64, (3, 2*N_events))

    for i in 1:N_events
        s = phi_summarize(lhe[i])
        res[:,2*i-1] = s[:,1]
        res[:,2*i] = s[:,2]
    end

    return res
end

function get_rts(rundir)
    n_run = rundir[end-2:end-1]
    if n_run == "00"
        n_run = rundir[end-3:end-1]
    end
    f = open(rundir * "run_$(n_run)_tag_1_banner.txt")
    l = ""
    for line in eachline(f)
        if occursin("ebeam1",line)
            l *= line
            break
        end
    end
    close(f)
    return 2*parse(Float64, strip(split(l,'=')[1]))
end

function get_mphi(rundir)
    if isfile(rundir * "params.dat")
        f = open(rundir * "params.dat")
        l = ""
        for line in eachline(f)
            if occursin("msca",line)
                l *= line
                break
            end
        end
        close(f)
        return parse(Float64, strip(split(l,"=")[2]))
    else
        return 0.0
    end
end

get_events(rundir) = parse_lhe(rundir * LHE_FILENAME)

function get_xsec(rundir)
    n_run = rundir[end-2:end-1]
    
    if n_run == "00"
        n_run = rundir[end-3:end-1]
        xsec_fname = rundir[1:end-15] * XSEC_RESULTS_FILENAME
    else
        xsec_fname = rundir[1:end-14] * XSEC_RESULTS_FILENAME
    end
    
    n_run = parse(Int64, n_run)
    return readdlm(xsec_fname, skipstart=1)[n_run,3]
end

function count_runs(mg_proc_dir)
    n_run = 0
    for ff in readdir(mg_proc_dir * "Events/")
        if isdir(mg_proc_dir * "Events/" * ff)
            n_run += 1
        end
    end
    return n_run
end

function get_all(fn, mg_proc_dir)
    n_run = count_runs(mg_proc_dir)

    placehold = Vector{Any}()

    for i in 1:n_run
        rundir = dir_runx(i)
        push!(placehold, fn(mg_proc_dir * rundir))
    end

    results = Vector{typeof(placehold[1])}(placehold)

    return results
end

function write_ll_summaries(dir)
    N_runs = count_runs(dir)

    for i in 1:N_runs
        rundir = dir*dir_runx(i)
        events_run = get_events(rundir)
        summary_run = dilepton_summarize(events_run)
        summary_df = DataFrame("m_ll" => summary_run[1,:], 
                               "pt_ll" => summary_run[2,:],
                               "cos_theta_ll" => summary_run[3,:],
                               "deltaR_ll" => summary_run[4,:],
                               "MT2" => summary_run[5,:])
        CSV.write(rundir * "dilepton_summary.csv", summary_df)
        
        println("$i of $(N_runs) complete.")
    end
end

function write_φφ_summaries(dir)
    N_runs = count_runs(dir)

    for i in 1:N_runs
        rundir = dir*dir_runx(i)
        events_run = get_events(rundir)
        summary_run = phi_summarize(events_run)
        summary_df = DataFrame("beta" => summary_run[1,:], 
                               "betagamma" => summary_run[2,:],
                               "eta" => summary_run[3,:])
        CSV.write(rundir * "diphi_summary.csv", summary_df)
        
        println("$i of $(N_runs) complete.")
    end
end

function write_ll_info(dir)
    rts = get_all(get_rts, dir)
    xsec = get_all(get_xsec, dir)
    df = DataFrame("rts" => rts, "xsec" => xsec)
    CSV.write(dir*"run_info.csv", df)

    write_ll_summaries(dir)
end

function write_ll_info(dir_χχ, dir_φφ)
    rts = get_all(get_rts, dir_χχ)
    xsec = get_all(get_xsec, dir_φφ)
    mphi = get_all(get_mphi, dir_χχ)
    df = DataFrame("rts" => rts, "mphi" => mphi, "xsec" => xsec)
    CSV.write(dir_χχ*"run_info.csv", df)

    write_ll_summaries(dir_χχ)
end

function write_φφ_info(dir)
    rts = get_all(get_rts, dir)
    xsec = get_all(get_xsec, dir)
    mphi = get_all(get_mphi, dir)
    df = DataFrame("rts" => rts, "mphi" => mphi, "xsec" => xsec)
    CSV.write(dir*"run_info.csv", df)

    write_φφ_summaries(dir)
end

function main()
    println("Background:")
    write_ll_info(MG_RESULTS_DIR*DIR_llνν)

    for f in FOLDERS[1:2]
        println("Prompt mumu:")
        write_ll_info(f*DIR_llχχ, f*DIR_φφ)

        println("Prompt VBF:")
        write_ll_info(f*DIR_VVχχ, f*DIR_VVφφ)
    end

    for f in FOLDERS[3:4]
        println("LLP mumu:")
        write_φφ_info(f*DIR_φφ)

        println("LLP VBF:")
        write_φφ_info(f*DIR_VVφφ)
    end
end

main()