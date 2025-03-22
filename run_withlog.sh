#!/bin/bash

# define a arry 
#puzzles=(
#"PZ1" "PZ2" "PZ3" "PZ4" "PZ5" "PZ6" "PZ7" "PZ8" "PZ9" "PZ10" "PZ12" "PZ13" "PZ14Bound" "PZ14Free" "PZ15" "PZ10tBox" "PZ10tRNA" "PZ11NMR" "PZ17" "PZ19" "PZ20" "PZ21"
#"PZ11" "PZ16a" "PZ16b" "PZ30"
#"PZ23" "PZ32"
#"PZ18" "PZ24" "PZ31" "PZ39"
#"PZ22" "PZ22Dimer" "PZ34" "PZ35" "PZ36"
#"PZ25" "PZ26" "PZ26tBox" "PZ26tRNA" "PZ27" "PZ27tBox" "PZ27tRNA" "PZ28" "PZ28tBox" "PZ28tRNA"
#"PZ29" "PZ33" "PZ37" "PZ38"
#)


#~ puzzles=(
#~ "PZ11" "PZ16a" "PZ16b" "PZ30"
#~ "PZ23" "PZ32"
#~ "PZ18" "PZ24" "PZ31" "PZ39"
#~ "PZ22" "PZ22Dimer" "PZ34" "PZ35" "PZ36"
#~ "PZ25" "PZ26" "PZ26tBox" "PZ26tRNA" "PZ27" "PZ27tBox" "PZ27tRNA" "PZ28" "PZ28tBox" "PZ28tRNA"
#~ "PZ29" "PZ33" "PZ37" "PZ38"
#~ )

#~ #"PZ36" "PZ35" "PZ33"  "PZ16a"
#~ puzzles=(
#~ "PZ29" "PZ30"
 #~ )
 
#puzzles=(
#"R0251" "R0254" "R0285" "R1248"
#)

#puzzles=(
#"R0281" "R0290" "R1203" "R1211" "R1221s2" "R1221s3" "R1224s2" "R1224s3" "R1241" "R1255" "R1256" "R1261" "R1262" "R1263" "R1264" "R1271" "R1289" "R1291" "R1209" "R1296" "R0250" "R0251" "R0252" "R0253" "R0254" "R0283" "R0285" "R1205" "R1212" "R1242" "R1248" "R1286" "R1288" "R1293"
#)
puzzles=("R1248" "R1288") 

DATE=20241123_casp16again


# touch a file of log
mkdir -p log_all/logout_${DATE}

# run script
for pz in "${puzzles[@]}"; do
    echo "Running for $pz..."
    
    
    # stdout and stderr all ---> log_out
    if ! python rnapuzzles_assess.py $pz - > "log_all/logout_${DATE}/${pz}_output.log" 2>&1; then
        echo "Error occurred while running for $pz. Check log_all/logout_${DATE}/${pz}_output.log for details."
    fi
done
