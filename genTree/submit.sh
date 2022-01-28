#!/bin/bash
date

if [ $# -ne 2 ]; then
     echo "Please input one arguement StartJobID, StopJobID ( (STOPID-STARID)should be <1000) !"
	 exit 1
fi

outdir=/star/data01/pwg/geurts/cc2ee/BES1_62GeV/out4ccbar

#rm -rf script/*
#rm -rf log/*

mkdir -p $outdir/output
mkdir -p $outdir/script
mkdir -p $outdir/log

ln -sf $outdir/output .
ln -sf $outdir/script .
ln -sf $outdir/log .

cp runMc.con runMcAll.con

jobStartID=$1
jobStopID=$2

while [ $jobStartID -le $jobStopID ]
do
  echo $jobStartID
  cp ./run.csh ./run_$jobStartID.csh

  echo "root -l -b <<EOF ">>run_$jobStartID.csh
  echo ".O2" >> run_$jobStartID.csh 
  echo ".L genPythia_C.so">> run_$jobStartID.csh
  echo -n "genPythia(">>run_$jobStartID.csh
  echo -n $jobStartID>>run_$jobStartID.csh
  echo -n ",5000000,">>run_$jobStartID.csh
  echo -n $RANDOM>>run_$jobStartID.csh
  echo -n ",0">>run_$jobStartID.csh
  echo ")">>run_$jobStartID.csh
  echo ".q">>run_$jobStartID.csh
  echo "EOF">>run_$jobStartID.csh

  #qsub -o log/run_$jobStartID.log -e log/run_$jobStartID.err run_$jobStartID.csh

  mv run_$jobStartID.csh script/

  echo "Executable     = /bin/sh">>runMcAll.con
  echo "Arguments      =" '"'"-c 'exec ./script/run_$jobStartID.csh'"'"'>>runMcAll.con
  echo "Output         = log/run_${jobStartID}.out" >>runMcAll.con
  echo "Error          = log/run_${jobStartID}.err" >>runMcAll.con
  echo "Log            = log/run_${jobStartID}.log" >>runMcAll.con
  echo "Queue"         >> runMcAll.con
  echo  "     " >>runMcAll.con

  let "jobStartID=jobStartID+1"

done

condor_submit runMcAll.con
