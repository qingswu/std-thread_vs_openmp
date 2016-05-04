#!/bin/bash

SIZES=("10000000" "100000000" "1000000000")
PROCS=("2" "4" "8" "16" "32")
for SIZE in "${SIZES[@]}"; do
    for THCT in "${PROCS[@]}"; do
	echo "###############################################################"
        ./bin/STREAM $THCT $SIZE
    done
done
