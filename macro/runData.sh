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

nEvents=1000
inputList=$1
nSkip=0

echo running: runData.sh $*
echo Running particle reconstruction with Fun4All_twoTrackReco.C
root.exe -q -b Fun4All_twoTrackReco.C\(${nEvents},\"${inputList}\",${nSkip}\)

if [[ "${useScratch}" == true ]]; then
  echo copying results back to $this_dir
  if [[ -d "output" ]]; then
    rsync -av "output" "${this_dir}/"
  fi
fi

echo Script done

# Note, to split a single run list into multiple run lists do gsplit -l 100 -d --additional-suffix=.txt run79516.list run79516_ (remove the first g on linux)
