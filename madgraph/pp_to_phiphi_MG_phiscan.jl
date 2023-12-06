# Parameters for Madgraph
mphis = Vector(100:50:2000) # GeV
mchi = 1e-6 # GeV
lambda = 1e-7
half_rts = 6500 # GeV

MADGRAPH_EXE = # insert path to Madgraph executable
MADGRAPH_OUTPUT = # insert location where you want to save madgraph outputs

DIR_φφ = "pp_to_phiphi_2"

script_file = MADGRAPH_OUTPUT*DIR_φφ*"_mgscript.txt"

function φφ_runcard(ebeam)
    φφ_stuff = """launch $DIR_φφ
    set mchi $mchi
    set msca scan:$(mphis)
    set lambda1 $lambda
    set lambda2 $lambda
    set lambda3 $lambda
    set nevents 100000
    set lpp1 1
    set lpp2 1
    set ebeam1 $ebeam
    set ebeam2 $ebeam
    0
    """
    return φφ_stuff
end

script_LFDM = """import model LFDM_UFO
generate p p > sca scabar QED=3
output $DIR_φφ
$(φφ_runcard(half_rts))launch $DIR_φφ -i
print_results --path=./$DIR_φφ/xsec.txt --format=short
exit
exit"""

if isdir(MADGRAPH_OUTPUT .* DIR_φφ)
    rm(MADGRAPH_OUTPUT .* DIR_φφ, recursive=true)
end

write(script_file, script_LFDM)
run(Cmd(`$MADGRAPH_EXE $script_file`, dir=MADGRAPH_OUTPUT))


cp(MADGRAPH_OUTPUT*DIR_φφ, (@__DIR__) * "/" * DIR_φφ, force=true)
