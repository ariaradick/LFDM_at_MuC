# Parameters for Madgraph
mphis = [Vector(range(1e3, 1.4e3, 50)), Vector(range(1e3, 4.9e3, 50))] # GeV
mchi = 1e-6 # GeV
lambda = 1e-7
half_rts = [1.5e3, 5e3] # GeV
RERUN_SM = false

MADGRAPH_EXE = # insert path to Madgraph executable
MADGRAPH_OUTPUT = # insert location where you want to save madgraph outputs

DIR_llνν = "mumu_to_llvv"
DIR_φφ = "mumu_to_phiphi"
DIR_llχχ = "mumu_to_llchichi"
DIR_VVφφ = "VV_to_phiphi"
DIR_VVllχχ = "VV_to_llchichi"

sm_script_file = MADGRAPH_OUTPUT*DIR_llνν*"_mgscript.txt"
lfdm_script_file = MADGRAPH_OUTPUT*DIR_llχχ*"_mgscript.txt"
VV_lfdm_script_file = MADGRAPH_OUTPUT*DIR_VVllχχ*"_mgscript.txt"

function sm_runcard(ebeams)
    SM_stuff = ""
    for ebeam in ebeams
        SM_stuff *= """launch $DIR_llνν
        set nevents 100000
        set lpp1 0
        set lpp2 0
        set ebeam1 $ebeam
        set ebeam2 $ebeam
        set sde_strategy 1
        0
        """
    end
    return SM_stuff
end

function φφ_runcard(ebeams)
    φφ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        φφ_stuff *= """launch $DIR_φφ
        set mchi $mchi
        set msca scan:$(mphis[i])
        set lambda1 $lambda
        set lambda2 $lambda
        set lambda3 $lambda
        set nevents 100000
        set lpp1 0
        set lpp2 0
        set ebeam1 $ebeam
        set ebeam2 $ebeam
        0
        """
    end
    return φφ_stuff
end

function χχ_runcard(ebeams)
    χχ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        χχ_stuff *= """launch $DIR_llχχ
        set mchi $mchi
        set msca scan:$(mphis[i])
        set lambda1 $lambda
        set lambda2 $lambda
        set lambda3 $lambda
        set nevents 100000
        set lpp1 0
        set lpp2 0
        set ebeam1 $ebeam
        set ebeam2 $ebeam
        set sde_strategy 1
        0
        """
    end
    return χχ_stuff
end

function VVφφ_runcard(ebeams)
    φφ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        φφ_stuff *= """launch $DIR_VVφφ
        set mchi $mchi
        set msca scan:$(mphis[i])
        set lambda1 $lambda
        set lambda2 $lambda
        set lambda3 $lambda
        set nevents 100000
        set lpp1 -4
        set lpp2 4
        set ebeam1 $ebeam
        set ebeam2 $ebeam
        set dsqrt_q2fact1 $ebeam
        set dsqrt_q2fact2 $ebeam
        0
        """
    end
    return φφ_stuff
end

function VVllχχ_runcard(ebeams)
    χχ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        χχ_stuff *= """launch $DIR_VVllχχ
        set mchi $mchi
        set msca scan:$(mphis[i])
        set lambda1 $lambda
        set lambda2 $lambda
        set lambda3 $lambda
        set nevents 100000
        set lpp1 -4
        set lpp2 4
        set ebeam1 $ebeam
        set ebeam2 $ebeam
        set dsqrt_q2fact1 $ebeam
        set dsqrt_q2fact2 $ebeam
        set sde_strategy 1
        0
        """
    end
    return χχ_stuff
end

script_SM = """generate mu+ mu- > l+ l- vl vl~
output $DIR_llνν
$(sm_runcard(half_rts))launch $DIR_llνν -i
print_results --path=./$DIR_llνν/xsec.txt --format=short
exit
exit"""

script_LFDM = """import model LFDM_UFO
generate mu+ mu- > sca scabar QED=3
output $DIR_φφ
$(φφ_runcard(half_rts))launch $DIR_φφ -i
print_results --path=./$DIR_φφ/xsec.txt --format=short
exit
generate mu+ mu- > sca scabar > l+ l- chi chi QED=3
output $DIR_llχχ
$(χχ_runcard(half_rts))launch $DIR_llχχ -i
print_results --path=./$DIR_llχχ/xsec.txt --format=short
exit
exit"""

script_LFDM_VV = """import model LFDM_UFO
set group_subprocesses false
define vxp = w+ z a
define vxm = w- z a
generate vxp vxm > sca scabar QED=3
output $DIR_VVφφ
$(VVφφ_runcard(half_rts))launch $DIR_VVφφ -i
print_results --path=./$DIR_VVφφ/xsec.txt --format=short
exit
generate vxp vxm > sca scabar > l+ l- chi chi QED=3
output $DIR_VVllχχ
$(VVllχχ_runcard(half_rts))launch $DIR_VVllχχ -i
print_results --path=./$DIR_VVllχχ/xsec.txt --format=short
exit
exit"""

for ff in (sm_script_file, lfdm_script_file)
    if isfile(ff)
        rm(ff)
    end
end

for d in (MADGRAPH_OUTPUT .* [DIR_φφ, DIR_llχχ, DIR_VVφφ, DIR_VVllχχ])
    if isdir(d)
        rm(d, recursive=true)
    end
end

if !(isdir(MADGRAPH_OUTPUT * DIR_llνν)) || RERUN_SM
    write(sm_script_file, script_SM)
    run(Cmd(`$MADGRAPH_EXE $sm_script_file`, dir=MADGRAPH_OUTPUT))
end

write(lfdm_script_file, script_LFDM)
run(Cmd(`$MADGRAPH_EXE $lfdm_script_file`, dir=MADGRAPH_OUTPUT))

write(VV_lfdm_script_file, script_LFDM_VV)
run(Cmd(`$MADGRAPH_EXE $VV_lfdm_script_file`, dir=MADGRAPH_OUTPUT))

for d in (DIR_llνν, DIR_φφ, DIR_llχχ, DIR_VVφφ, DIR_VVllχχ)
    cp(MADGRAPH_OUTPUT*d, (@__DIR__) * "/phiscan/" * d, force=true)
end
