#!/bin/sh 
# $1 - sample

#/eos/cms/store/group/phys_tau/TauML/prod_2018_v2/full_tuples/DYJetsToLL_M-50
FILENAME=${1##*/}
FILENAME=${FILENAME%.*}
cat > jobs/runTuple_$FILENAME.submit <<EOF
+RequestRuntime=10000

RequestMemory = 2000

executable = jobs/runTuple_$FILENAME.sh
arguments = $FILENAME output

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "EL7"
when_to_transfer_output = ON_EXIT

output              = jobs/$FILENAME.out
error               = jobs/$FILENAME.error
log                 = jobs/$FILENAME.log

queue

EOF

cat > jobs/runTuple_$FILENAME.sh <<EOF1
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd ${CMSSW_BASE}/src
cmsenv
cd -
runTuple $FILENAME output
EOF1
chmod u+x jobs/runTuple_$FILENAME.sh
chmod u+x jobs/runTuple_$FILENAME.submit
condor_submit jobs/runTuple_$FILENAME.submit
