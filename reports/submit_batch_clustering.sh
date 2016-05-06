#!/bin/bash

qsub -cwd -b y -q bhklab -e sge_out -o sge_out -N "cluster" -t 1-4500 -tc 50 "module load R; Rscript batch.cluster.all.R"
qsub -cwd -b y -q bhklab -e sge_out -o sge_out -N "cluster" -t 1-4500 -tc 50 "module load R; Rscript batch.cluster.concordant.R"
