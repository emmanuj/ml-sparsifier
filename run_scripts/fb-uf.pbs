#!/bin/bash
#PBS -N fb-uf
#PBS -l select=1:ncpus=16:mem=60gb
#PBS -l walltime=7:00:00
#PBS -M emmanuj@g.clemson.edu
cd $PBS_O_WORKDIR
cd ..
module load gcc/4.8.1 python/3.4
./msparse -i ~/test_graphs/fb-uf.edges -n 35111 --weak --zero-index -o out/fb-uf_coarsest.txt -r 1 --level-span 3 --sparse-level 0 -s 0.3
./msparse -i ~/test_graphs/fb-uf.edges -n 35111 --weak --zero-index -o out/fb-uf_mid.txt -r 1 --level-span 3 --sparse-level 1 -s 0.3
./msparse -i ~/test_graphs/fb-uf.edges -n 35111 --weak --zero-index -o out/fb-uf_finest.txt -r 1 --level-span 3 --sparse-level 2 -s 0.3
python3 benchmark.py fb-uf /home/emmanuj/test_graphs/fb-uf.edges
