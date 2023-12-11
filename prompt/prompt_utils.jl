module Prompt

export load_run_info, load_summary, get_bkg_at_rts, dir_runx,
    significances, fcut, N_signal, N_background

using CSV
using DataFrames

RUN_RESULTS_FILENAME = "run_info.csv"
SUMMARY_FILENAME = "dilepton_summary.csv"

function dir_runx(x)
    if x < 10
        return "Events/run_0$x/"
    else
        return "Events/run_$x/"
    end
end

function load_summary(rundir)
    Matrix(CSV.read(rundir*SUMMARY_FILENAME, DataFrame))'
end

function load_run_info(mgdir)
    CSV.read(mgdir*RUN_RESULTS_FILENAME, DataFrame)
end

function get_bkg_at_rts(rts_tev, bkg_dir)
    rts_gev = Float64(rts_tev)*1e3

    run_info = load_run_info(bkg_dir)
    idxs = findall(x -> x == rts_gev, run_info[:,1])
    run_info = run_info[idxs,:]

    if length(idxs) > 1
        throw(ArgumentError("too many bkgs"))
    else
        bkg_summary = load_summary(bkg_dir*dir_runx(idxs[1]))
    end

    return run_info, bkg_summary
end

gtr_min(summ_vec, min) = summ_vec .> min

function fcut(summary, mll_min, ptll_min, mt2_min)
    N_mll = length(mll_min)
    N_ptll = length(ptll_min)
    N_mt2 = length(mt2_min)
    Nevent = length(summary[1,:])

    res = zeros(Float64, (N_mll, N_ptll, N_mt2))

    mll_test = gtr_min.((summary[1,:],), mll_min)
    ptll_test = gtr_min.((summary[2,:],), ptll_min)
    mt2_test = gtr_min.((summary[5,:],), mt2_min)

    @Threads.threads for i in 1:N_mll
        @Threads.threads for j in 1:N_ptll
            @Threads.threads for k in 1:N_mt2
                res[i,j,k] = count(mll_test[i] .&& ptll_test[j] .&& mt2_test[k])
            end
        end
    end

    return res ./ Nevent
end

fcut(summary, mins) = fcut(summary, mins...)[1]

function N_signal(run, lumi, mins)
    N_mumu = run.mumu_xsec * lumi * fcut(run.mumu_summary, mins)
    N_vv = run.vv_xsec * lumi * fcut(run.vv_summary, mins)

    return N_mumu+N_vv
end

function N_background(run, lumi, mins)
    run.bkg_xsec * lumi * fcut(run.bkg_summary, mins)
end

function significances(run, lumi, mll_min, ptll_min, mt2_min)
    mumu_frac = fcut(run.mumu_summary, mll_min, ptll_min, mt2_min)
    vv_frac = fcut(run.vv_summary, mll_min, ptll_min, mt2_min)
    bkg_frac = fcut(run.bkg_summary, mll_min, ptll_min, mt2_min)

    N_mumu = run.mumu_xsec .* lumi .* mumu_frac
    N_vv = run.vv_xsec .* lumi .* vv_frac
    N_bkg_true = run.bkg_xsec .* lumi .* bkg_frac
    N_bkg = N_bkg_true .+ ((N_bkg_true .<= 2.0) .* (2.0 .- N_bkg_true))

    N_S = @. (N_mumu + N_vv)

    @. N_S *= (N_S >= 20)

    return N_S ./ sqrt.(N_bkg)
end

end