#!/bin/bash

PREFIX="simple_rcb"
ptRange=$(awk 'BEGIN{for(i=5;i<=60;i+=0.1)print i}')
threshRange=$(awk 'BEGIN{for(j=0.1;j<=0.3;j+=0.1)print j}')

echo "Simulating $PARTICLE events."
# Require that user specify particle type. 
if [ $# -eq 0 ]; then 
    echo "Error: didn't specify EMinus, Gamma, or Pi0. Exiting."
    exit
elif [ $# -eq 1 ] && [ "$1" == "test" ]; then 
    echo "Starting a test run . . . "
    cd ~/macros/g4simulations;
    root -b -q "Fun4All_G4_sPHENIX.C(1, 10.0, 0.2, \"EMinus\", 0)"
    exit
fi

# Store desired simulated particle type.
PARTICLE="$1"
 
# Run single-particle events for genPT from 5 - 60GeV. 
cd ~/macros/g4simulations;
EVENT_NUM=0
for seedThresh in $threshRange; do

    EVENT_NUM=0

    for genPT in $ptRange; do
        root -b -q "Fun4All_G4_sPHENIX.C(1, $genPT, $seedThresh, \"${PARTICLE}\", $EVENT_NUM)"
        echo "Finished event with genPT = $genPT . . . "
        EVENT_NUM=$((EVENT_NUM+1))
    done

    # Combine root files 
    # and place the leftover individual files in their own folder.
    cd ~/bmckinz/RawClusterBuilderIA/rootFiles
    TARGET=${PREFIX}_${seedThresh}_${PARTICLE}.root
    if [ -f $TARGET ]; then
        echo "Removing existing $TARGET . . . "
        rm $TARGET
    fi
    hadd ${TARGET} *${PARTICLE}_*
    mv *${PARTICLE}_* $PARTICLE/
    cd -

done

cd ~/bmckinz/RawClusterBuilderIA/rootFiles
TARGET=${PREFIX}${PARTICLE}.root
if [ -f $TARGET ]; then
    echo "Removing existing $TARGET . . . "
    rm $TARGET
fi
hadd ${TARGET} *_${PARTICLE}.root
mv *_${PARTICLE}.root $PARTICLE/

cd -;

