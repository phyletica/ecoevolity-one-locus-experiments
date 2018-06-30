#! /bin/bash

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR
fi

nsites="02000"

nprocs=8
nprior=500000
batch_size=6250
nsums=100000
seed=1384268

output_dir="../prior-${nsites}"
if [ ! -d "$output_dir" ]
then
    mkdir -p $output_dir
fi

dmc.py --np $nprocs -r 1 -o ../configs/dpp-msbayes-${nsites}.cfg -p ../configs/dpp-msbayes-${nsites}.cfg -n $nprior --num-posterior-samples $batch_size --prior-batch-size $batch_size --num-standardizing-samples $nsums --output-dir $output_dir --seed $seed --generate-samples-only 1>generate-prior-samples-${nsites}.sh.out 2>&1
