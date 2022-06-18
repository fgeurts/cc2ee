#!/bin/bash

outDir=/star/data01/pwg/geurts/cc2ee/BES1_27GeV/out4cchist

mkdir -p ${outDir}/csh
mkdir -p ${outDir}/report
mkdir -p ${outDir}/list
mkdir -p ${outDir}/log
mkdir -p ${outDir}/err
mkdir -p ${outDir}/output

ln -sf ${outDir} .
ln -sf ${outDir}/csh .
ln -sf ${outDir}/report .
ln -sf ${outDir}/list .
ln -sf ${outDir}/log .
ln -sf ${outDir}/err .
ln -sf ${outDir}/output .

dir=$(pwd)
echo "current directory: $dir"

ioutDirCent=${outDir}/output
mkdir -p ${ioutDirCent}
star-submit-template -template submit4cchists.xml -entities logpath=${outDir},outpath=${ioutDirCent}
