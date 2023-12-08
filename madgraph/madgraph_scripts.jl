# Run this to generate the madgraph data you'll need to run the other parts of
# this code. Please make sure to insert the path to your madgraph executable,
# and the output folder you wish madgraph to use (cannot contain spaces).
# This script will automatically move the madgraph output directories to the 
# location of this file.

MADGRAPH_EXE = "/home/aradick/Downloads/Madgraph/MG5_aMC_v3_5_1/bin/mg5_aMC" # insert path to Madgraph executable
MADGRAPH_OUTPUT = "/home/aradick/Downloads/Madgraph/test/" # insert location where you want to save madgraph outputs

# Parameters for Madgraph
mphis_prompt = [[1e3, 1.3e3], [1e3, 4.5e3]] # GeV
mphis_prompt_scan = [Vector(range(1e3, 1.4e3, 50)), Vector(range(1e3, 4.9e3, 50))] # GeV

mphis_LLP = [[500, 1e3, 1.3e3], [500, 1e3, 4.5e3]] # GeV
mphis_LLP_scan = [Vector(10 .^ (range(2, log10(1.5e3), 51)[1:end-1])), 
                  Vector(10 .^ (range(2, log10(5e3), 51)[1:end-1]))] # GeV

Mphis = [mphis_prompt, mphis_prompt_scan, mphis_LLP, mphis_LLP_scan]

mchi = 1e-6 # GeV
lambda = 1e-7
half_rts = [1.5e3, 5e3] # GeV

DIR_llνν = "mumu_to_llvv" # output dir for background, mumu to llνν

DIR_φφ = "mumu_to_phiphi" # output dir for mumu to φφ, necessary to get the cross-section
DIR_llχχ = "mumu_to_llchichi" # output dir for mumu to llχχ

DIR_VVφφ = "VV_to_phiphi" # output dir for VBF to φφ, necessary to get the cross-section
DIR_VVllχχ = "VV_to_llchichi" # output dir for VBF to llχχ

DIR_AZφφ = "AZ_to_phiphi"
DIR_WWφφ = "WW_to_phiphi"

FOLDERS = ["prompt/", "prompt_scan/", "LLP/", "LLP_scan/", "xsec/"]
do_χχ = [true, true, false, false]

# names and locations of the madgraph scripts that will be generated and run
sm_script_file = MADGRAPH_OUTPUT*DIR_llνν*"_mgscript.txt"
lfdm_script_file = DIR_llχχ*"_mgscript.txt"
VV_lfdm_script_file = DIR_VVllχχ*"_mgscript.txt"

AZ_script_file = MADGRAPH_OUTPUT*"xsec/"*DIR_AZφφ*"_mgscript.txt"
WW_script_file = MADGRAPH_OUTPUT*"xsec/"*DIR_WWφφ*"_mgscript.txt"

# functions to generate the part of the MG script that changes the run card
# and param card, necessary to scan over rts values.
function sm_script(ebeams)
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

    script = """generate mu+ mu- > l+ l- vl vl~
    output $DIR_llνν
    $(SM_stuff)launch $DIR_llνν -i
    print_results --path=./$DIR_llνν/xsec.txt --format=short
    exit
    exit"""

    return script
end

function μμ_script(ebeams, mass_phis; χχ=true)
    φφ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        φφ_stuff *= """launch $DIR_φφ
        set mchi $mchi
        set msca scan:$(mass_phis[i])
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
    
    script = """import model LFDM_UFO
    generate mu+ mu- > sca scabar QED=3
    output $DIR_φφ
    $(φφ_stuff)launch $DIR_φφ -i
    print_results --path=./$DIR_φφ/xsec.txt --format=short
    exit
    """

    if χχ
        χχ_stuff = ""
        for (i,ebeam) in enumerate(ebeams)
            χχ_stuff *= """launch $DIR_llχχ
            set mchi $mchi
            set msca scan:$(mass_phis[i])
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

        script *= """generate mu+ mu- > sca scabar > l+ l- chi chi QED=3
        output $DIR_llχχ
        $(χχ_stuff)launch $DIR_llχχ -i
        print_results --path=./$DIR_llχχ/xsec.txt --format=short
        exit
        exit"""
    else
        script *= "exit"
    end

    return script
end

function VV_script(ebeams, mass_phis; χχ=true)
    φφ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        φφ_stuff *= """launch $DIR_VVφφ
        set mchi $mchi
        set msca scan:$(mass_phis[i])
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

    script = """import model LFDM_UFO
    set group_subprocesses false
    define vxp = w+ z a
    define vxm = w- z a
    generate vxp vxm > sca scabar QED=3
    output $DIR_VVφφ
    $(φφ_stuff)launch $DIR_VVφφ -i
    print_results --path=./$DIR_VVφφ/xsec.txt --format=short
    exit
    """

    if χχ
        χχ_stuff = ""
        for (i,ebeam) in enumerate(ebeams)
            χχ_stuff *= """launch $DIR_VVllχχ
            set mchi $mchi
            set msca scan:$(mass_phis[i])
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

        script *= """generate vxp vxm > sca scabar > l+ l- chi chi QED=3
        output $DIR_VVllχχ
        $(χχ_stuff)launch $DIR_VVllχχ -i
        print_results --path=./$DIR_VVllχχ/xsec.txt --format=short
        exit
        exit"""
    else
        script *= "exit"
    end

    return script
end

function VVφφ_runcard(dir, ebeams)
    φφ_stuff = ""
    for (i,ebeam) in enumerate(ebeams)
        φφ_stuff *= """launch $dir
        set mchi $mchi
        set msca scan:$(mphis_LLP_scan[i])
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

# write the scripts to files, then run with madgraph:

write(sm_script_file, sm_script(half_rts))
run(Cmd(`$MADGRAPH_EXE $sm_script_file`, dir=MADGRAPH_OUTPUT))

for (i, d) in enumerate(FOLDERS)
    if i == 5
        break
    end

    ff = MADGRAPH_OUTPUT*d

    mkdir(ff)

    write(ff*lfdm_script_file, μμ_script(half_rts, Mphis[i]; χχ=do_χχ[i]))
    run(Cmd(`$MADGRAPH_EXE $lfdm_script_file`, dir=ff))

    write(ff*VV_lfdm_script_file, VV_script(half_rts, Mphis[i]; χχ=do_χχ[i]))
    run(Cmd(`$MADGRAPH_EXE $VV_lfdm_script_file`, dir=ff))
end

mkdir(MADGRAPH_OUTPUT*"xsec/")

write(AZ_script_file, script_AZ)
run(Cmd(`$MADGRAPH_EXE $AZ_script_file`, dir=MADGRAPH_OUTPUT*"xsec/"))

write(WW_script_file, script_WW)
run(Cmd(`$MADGRAPH_EXE $WW_script_file`, dir=MADGRAPH_OUTPUT*"xsec/"))

# copy output dirs to the location of this folder
for d in FOLDERS
    cp(MADGRAPH_OUTPUT*d, (@__DIR__) * "/" * d, force=true)
end