#!/bin/bash

NTASKS=("1" "2" "4" "8" "16" "32")
PROCS=("2" "4" "8" "16" "32")
for THCT in "${PROCS[@]}"; do
    for NTASK in "${NTASKS[@]}"; do
        echo "##################################################"
        #./bin/streamLFTPntasks $THCT 100000000 $(($THCT*$NTASK))
        ./bin/streamTP11ntasks $THCT 100000000 $(($THCT*$NTASK))
    done
done
