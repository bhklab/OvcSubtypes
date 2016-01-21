#!/bin/bash

qsub -cwd -b y -q bhklab -e sge_out -o sge_out -N "cluster" -t 2016-4800 -tc 15 "module load R; Rscript batch.cluster.all.R"
