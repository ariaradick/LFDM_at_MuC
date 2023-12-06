using Pkg
Pkg.activate("../one_flav_env")

using LHEF
using LorentzVectorHEP
using DataFrames
using CSV
using StaticArrays
using LinearAlgebra
using DelimitedFiles: readdlm
using Optim

# MG_RESULTS_DIR = "../madgraph/phiscan/"
MG_RESULTS_DIR = "../madgraph/"

# make sure these match what the script has output
DIR_llνν = MG_RESULTS_DIR * "mumu_to_llvv/"
DIR_φφ = MG_RESULTS_DIR * "mumu_to_phiphi/"
DIR_llχχ = MG_RESULTS_DIR * "mumu_to_llchichi/"
DIR_VVφφ = MG_RESULTS_DIR * "VV_to_phiphi/"
DIR_VVχχ = MG_RESULTS_DIR * "VV_to_llchichi/"

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

function trans_mass(pT1, pT2, m1, m2)
    ET1 = sqrt((pT1 ⋅ pT1) + m1^2)
    ET2 = sqrt((pT2 ⋅ pT2) + m2^2)
    return m1^2 + m2^2 + 2*(ET1*ET2 - (pT1 ⋅ pT2))
end

function MT2(particle_1, particle_2)
    p1 = lhe_v4(particle_1)
    p2 = lhe_v4(particle_2)

    m1 = particle_1.m
    m2 = particle_2.m

    pT1 = lhe_pt(particle_1)
    pT2 = lhe_pt(particle_2)
    pt_miss = -(pT1 .+ pT2)
    init_x = MArray{Tuple{2}}(pt_miss) ./ 2

    function mt_to_min(x)
        ET1 = sqrt((pT1 ⋅ pT1) + m1^2)
        ET2 = sqrt((pT2 ⋅ pT2) + m2^2)
        return max(trans_mass(pT1, x, m1, 0), trans_mass(pT2, pt_miss-x, m2, 0))
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

function runs_summary(bkg_dir, sig_dirs, φφ_dirs)
    rts_bkg = get_all(get_rts, bkg_dir)
    sig_mumu = get_all.((get_rts, get_mphi), sig_dirs[1])
    sig_vv = get_all.((get_rts, get_mphi), sig_dirs[2])

    xsecs_bkg = get_all(get_xsec, bkg_dir)
    xsecs_mumu = get_all(get_xsec, φφ_dirs[1])
    xsecs_vv = get_all(get_xsec, φφ_dirs[2])

    if !(sig_mumu == sig_vv)
        throw(ArgumentError("mumu and vbf don't align!"))
    end

    x_bkg = DataFrame("rts" => rts_bkg,
                       "xsec" => xsecs_bkg)

    x_mumu = DataFrame("rts" => sig_mumu[1], 
                       "mphi" => sig_mumu[2],
                       "xsec" => xsecs_mumu)
    
    x_vv = DataFrame("rts" => sig_vv[1], 
                     "mphi" => sig_vv[2],
                     "xsec" => xsecs_vv)

    CSV.write(bkg_dir*"run_info.csv", x_bkg)
    CSV.write(sig_dirs[1]*"run_info.csv", x_mumu)
    CSV.write(sig_dirs[2]*"run_info.csv", x_vv)

    completion_text = ["bkg", "mumu", "vbf"]
    dirs = [bkg_dir, sig_dirs[1], sig_dirs[2]]
    N_runs = count_runs.(dirs)

    for j in 1:3
        for i in 1:N_runs[j]
            rundir = dirs[j]*dir_runx(i)
            events_run = get_events(rundir)
            summary_run = dilepton_summarize(events_run)
            summary_df = DataFrame("m_ll" => summary_run[1,:], "pt_ll" => summary_run[2,:],
                                    "cos_theta_ll" => summary_run[3,:],
                                    "deltaR_ll" => summary_run[4,:],
                                    "MT2" => summary_run[5,:])
            CSV.write(rundir * "dilepton_summary.csv", summary_df)

            events_run = nothing
            summary_run = nothing
            summary_df = nothing
            GC.gc()
            
            println("$(completion_text[j]): $i of $(N_runs[j]) complete.")
        end
    end
end

function main()
    runs_summary(DIR_llνν, (DIR_llχχ, DIR_VVχχ), (DIR_φφ, DIR_VVφφ))
end

main()