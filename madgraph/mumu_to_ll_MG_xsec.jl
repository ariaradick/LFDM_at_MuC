# Parameters for Madgraph
mphis = [Vector(10 .^ (range(2, log10(1.5e3), 51)[1:end-1])), 
         Vector(10 .^ (range(2, log10(5e3), 51)[1:end-1]))] # GeV
mchi = 1e-6 # GeV
lambda = 1e-7
half_rts = [1.5e3, 5e3] # GeV

MADGRAPH_EXE = # insert path to Madgraph executable
MADGRAPH_OUTPUT = # insert location where you want to save madgraph outputs

DIR_φφ = "mumu_to_phiphi_xsec"
DIR_AZφφ = "AZ_to_phiphi_xsec"
DIR_WWφφ = "WW_to_phiphi_xsec"

mumu_script_file = MADGRAPH_OUTPUT*DIR_φφ*"_mgscript.txt"
AZ_script_file = MADGRAPH_OUTPUT*DIR_AZφφ*"_mgscript.txt"
WW_script_file = MADGRAPH_OUTPUT*DIR_WWφφ*"_mgscript.txt"

function mumu_φφ_runcard(ebeams)
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

function VVφφ_runcard(dir, ebeams)
    φφ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        φφ_stuff *= """launch $dir
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

script_mumu = """import model LFDM_UFO
generate mu+ mu- > sca scabar QED=3
output $DIR_φφ
$(mumu_φφ_runcard(half_rts))launch $DIR_φφ -i
print_results --path=./$DIR_φφ/xsec.txt --format=short
exit
exit"""

script_AZ = """import model LFDM_UFO
set group_subprocesses false
generate a a > sca scabar QED=3
add process a z > sca scabar QED=3
add process z a > sca scabar QED=3
add process z z > sca scabar QED=3
output $DIR_AZφφ
$(VVφφ_runcard(DIR_AZφφ,half_rts))launch $DIR_AZφφ -i
print_results --path=./$DIR_AZφφ/xsec.txt --format=short
exit
exit"""

script_WW = """import model LFDM_UFO
set group_subprocesses false
generate w+ w- > sca scabar QED=3
output $DIR_WWφφ
$(VVφφ_runcard(DIR_WWφφ,half_rts))launch $DIR_WWφφ -i
print_results --path=./$DIR_WWφφ/xsec.txt --format=short
exit
exit"""

# for ff in (mumu_script_file, AZ_script_file, WW_script_file)
for ff in (AZ_script_file, WW_script_file)
    if isfile(ff)
        rm(ff)
    end
end

# for d in (MADGRAPH_OUTPUT .* [DIR_φφ, DIR_AZφφ, DIR_WWφφ])
for d in (MADGRAPH_OUTPUT .* [DIR_AZφφ, DIR_WWφφ])
    if isdir(d)
        rm(d, recursive=true)
    end
end

# write(mumu_script_file, script_mumu)
# run(Cmd(`$MADGRAPH_EXE $mumu_script_file`, dir=MADGRAPH_OUTPUT))

write(AZ_script_file, script_AZ)
run(Cmd(`$MADGRAPH_EXE $AZ_script_file`, dir=MADGRAPH_OUTPUT))

write(WW_script_file, script_WW)
run(Cmd(`$MADGRAPH_EXE $WW_script_file`, dir=MADGRAPH_OUTPUT))

for d in (DIR_φφ, DIR_AZφφ, DIR_WWφφ)
    cp(MADGRAPH_OUTPUT*d, (@__DIR__) * "/" * d, force=true)
end
