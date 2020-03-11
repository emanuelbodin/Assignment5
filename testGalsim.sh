#!/bin/bash

N=(500, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000)

Nfile=(00500 01000 02000 04000 06000 08000 10000 12000 14000 16000)

echo "Testing for 2 cores"

for i in {0..10}
do
    time ./galsim ${N[$i]} input_data/ellipse_N_${Nfile[$i]}.gal 200 1e-5 0.255 0 2
done

echo "Testing for 2 cores"

for i in {0..10}
do
    time ./galsim ${N[$i]} input_data/ellipse_N_${Nfile[$i]}.gal 200 1e-5 0.255 0 4
done

echo "Testing for 1 to 4 cores on N = 500"

for i in {0..4}
do
    time ./galsim ${N[0]} input_data/ellipse_N_${Nfile[0]}.gal 200 1e-5 0.255 0 $i
done

echo "Testing for 1 to 4 cores on N = 6000"

for i in {0..4}
do
    time ./galsim ${N[4]} input_data/ellipse_N_${Nfile[4]}.gal 200 1e-5 0.255 0 $i
done