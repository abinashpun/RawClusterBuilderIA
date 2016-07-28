#!/bin/bash

# Require that user specify particle type. 
if [ $# -eq 0 ]; then 
    echo "Error: didn't specify EMinus, Gamma, or Pi0. Exiting."
    exit
elif [ $# -eq 1 ] && [ "$1" == "test" ]; then 
    echo "Starting a test run . . . "
    cd ~/macros/g4simulations;
    root -b -q "Fun4All_G4_sPHENIX.C(1, 10.0, \"EMinus\", 0)"
    exit
fi

# Store desired simulated particle type.
PARTICLE="$1"
echo "Simulating $PARTICLE events."
#
# Run single-particle events for genPT from 5 - 60GeV. 
cd ~/macros/g4simulations;
EVENT_NUM=0
ptRange=$(awk 'BEGIN{for(i=5;i<=60;i+=1)print i}')
for genPT in $ptRange; do
    root -b -q "Fun4All_G4_sPHENIX.C(1, $genPT, \"${PARTICLE}\", $EVENT_NUM)"
    echo "Finished event with genPT = $genPT . . . "
    EVENT_NUM=$((EVENT_NUM+1))
done
cd -;

# Combine root files 
# and place the leftover individual files in their own folder.
cd ~/bmckinz/MyRawClusterBuilder/rootFiles
TARGET=rcb_${PARTICLE}.root
if [ -f $TARGET ]; then
    echo "Removing existing $TARGET . . . "
    rm $TARGET
fi
hadd rcb_${PARTICLE}.root *${PARTICLE}_*
mv *${PARTICLE}_* $PARTICLE/
cd -

