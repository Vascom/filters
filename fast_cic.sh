#!/bin/bash

g++ -O2 cic_comp.cpp -o cic3.bin
./cic3.bin 190 10 -100 -99 190_10_0_10000.txt 1 &
./cic3.bin 190 10 -98 -97 190_10_1_10000.txt 1 &
./cic3.bin 190 10 -96 -95 190_10_2_10000.txt 1 &
./cic3.bin 190 10 -94 -93 190_10_3_10000.txt 1 &
