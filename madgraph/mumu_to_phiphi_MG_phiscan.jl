# Parameters for Madgraph
# mphis = [Vector(10 .^ (range(2, log10(1.5e3), 51)[1:end-1])), 
#          Vector(10 .^ (range(2, log10(5e3), 51)[1:end-1]))] # GeV
mphis = [[500, 1e3, 1.3e3], [500, 1e3, 4.5e3]]
mchi = 1e-6 # GeV
lambda = 1e-7
half_rts = [1.5e3, 5e3] # GeV

MADGRAPH_EXE = # insert path to Madgraph executable
MADGRAPH_OUTPUT = # insert location where you want to save madgraph outputs

DIR_φφ = "mumu_to_phiphi_2"
DIR_VVφφ = "VV_to_phiphi_2"

lfdm_script_file = MADGRAPH_OUTPUT*DIR_φφ*"_mgscript.txt"
VV_lfdm_script_file = MADGRAPH_OUTPUT*DIR_VVφφ*"_mgscript.txt"

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

script_LFDM = """import model LFDM_UFO
generate mu+ mu- > sca scabar QED=3
output $DIR_φφ
$(φφ_runcard(half_rts))launch $DIR_φφ -i
print_results --path=./$DIR_φφ/xsec.txt --format=short
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
exit"""

for d in (MADGRAPH_OUTPUT .* [DIR_φφ, DIR_VVφφ])
# for d in (MADGRAPH_OUTPUT .* [DIR_VVφφ])
    if isdir(d)
        rm(d, recursive=true)
    end
end

write(lfdm_script_file, script_LFDM)
run(Cmd(`$MADGRAPH_EXE $lfdm_script_file`, dir=MADGRAPH_OUTPUT))

write(VV_lfdm_script_file, script_LFDM_VV)
run(Cmd(`$MADGRAPH_EXE $VV_lfdm_script_file`, dir=MADGRAPH_OUTPUT))

# for d in (DIR_φφ, DIR_VVφφ)
# for d in (DIR_VVφφ)
#     cp(MADGRAPH_OUTPUT*d, (@__DIR__) * "/" * d, force=true)
# end
