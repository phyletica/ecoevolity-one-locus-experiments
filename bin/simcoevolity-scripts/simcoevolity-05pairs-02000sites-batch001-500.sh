#! /bin/sh

if [ -n "$PBS_JOBNAME" ]
then
    source ${PBS_O_HOME}/.bash_profile
    cd $PBS_O_WORKDIR

    module load gcc/5.3.0
fi

locussize="2000"
simname="05pairs-0${locussize}sites"
cfgpath="../../configs/config-${simname}.yml"
outputdir="../../simulations/validation/${simname}/batch001"
rngseed=784052805
nreps=500

mkdir -p "$outputdir"

simcoevolity --seed="$rngseed" -n "$nreps" -l "$locussize" -o "$outputdir" "$cfgpath"
