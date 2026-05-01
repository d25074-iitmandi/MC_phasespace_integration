#!/bin/bash
export OMP_NUM_THREADS=1
make clean && make

rm -f results/scaling.txt

for t in 1 2 4 8
do
    export OMP_NUM_THREADS=$t
    ./mc
done
