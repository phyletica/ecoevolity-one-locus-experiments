#! /bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=120:00:00
#PBS -l mem=120gb
#PBS -j oe 
#PBS -l jobflags=ADVRES:jro0014_lab.56281

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR
fi

reps=500
nprocs=8
nprior=500000
batch_size=12500
nsums=100000
npost=2000
nquantiles=1000
sortindex=11
seed=37851841

output_dir="../results/dpp-msbayes"
if [ ! -d "$output_dir" ]
then
    mkdir -p $output_dir
fi

dmc.py --np $nprocs \
    -r $reps \
    -o ../configs/dpp-msbayes.cfg \
    -p ../prior/pymsbayes-results/pymsbayes-output/prior-stats-summaries \
    -n $nprior \
    --prior-batch-size $batch_size \
    --num-posterior-samples $npost \
    --num-standardizing-samples $nsums \
    -q $nquantiles \
    --sort-index $sortindex \
    --output-dir $output_dir \
    --seed $seed \
    --no-global-estimate \
    --compress \
    1>run-reject-dpp-msbayes.sh.out 2>&1
