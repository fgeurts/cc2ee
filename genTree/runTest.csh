#!/bin/bash
mkdir -p outtest
rm -rf outtest/test.log
rm -rf outtest/pythiaevent_test.root


#void genPythia(Int_t irun=1, Int_t Nevt=1000, Int_t iseed=789456, Bool_t debug=1){

root4star -b <<EOF >& outtest/test.log
.O2
.L genPythia_C.so
genPythia(888,50000,3128,1);
.q
EOF
