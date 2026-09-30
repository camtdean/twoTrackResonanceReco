#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new 

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}

# uncomment for local development
export SPHENIX=${HOME}/sPHENIX
export MYINSTALL=$SPHENIX/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH
export PATH="$HOME/.local/bin:$PATH"

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

useScratch=false

if [[ "${useScratch}" == true ]]; then
  this_script=$BASH_SOURCE
  this_script=`readlink -f $this_script`
  this_dir=`dirname $this_script`
  echo rsyncing from $this_dir
  echo running: $this_script $*
  
  if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
  then
    cd $_CONDOR_SCRATCH_DIR
    rsync -av $this_dir/* .
  else
    echo condor scratch NOT set
    exit -1
  fi
fi

nEvents=100
outDir=./
processID=$1

echo running: runSims.sh $*
echo D0 simulation with CA seeding
root.exe -q -b Fun4All_D0_sim.C\(${nEvents},\"${outDir}\",${processID},false\)
echo D0 simulation with polyseeding
root.exe -q -b Fun4All_D0_sim.C\(${nEvents},\"${outDir}\",${processID},true\)
echo Script done