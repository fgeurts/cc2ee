#!/bin/bash

#beamenergy=27,39,62
beamenergy=$1
outDir=/star/data01/pwg/geurts/cc2ee/BES1_${beamenergy}GeV/out4cchist

mkdir -p ${outDir}/csh
mkdir -p ${outDir}/report
mkdir -p ${outDir}/list
mkdir -p ${outDir}/log
mkdir -p ${outDir}/err
mkdir -p ${outDir}/output

dir=$(pwd)
echo "current directory: $dir"

ioutDirCent=${outDir}/output
mkdir -p ${ioutDirCent}
star-submit-template -template submit4cchists.xml -entities logpath=${outDir},outpath=${ioutDirCent},energy=${beamenergy}
