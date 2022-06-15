#!/bin/sh 
# $1 - full path input file

#/eos/cms/store/group/phys_tau/TauML/prod_2018_v2/full_tuples/DYJetsToLL_M-50
FILENAME=${1##*/}
FILENAME=${FILENAME%.*}
cat > jobs/runTuple_$FILENAME.submit <<EOF
+RequestRuntime = 1200

executable = jobs/runTuple_$FILENAME.sh

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "CentOS7"

output              = jobs/$FILENAME.out
error               = jobs/$FILENAME.error
log                 = jobs/$FILENAME.log

queue

EOF

cat > jobs/runTuple_$FILENAME.sh <<EOF1
#!/bin/sh
export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-pafccj 
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700  
cd ${CMSSW_BASE}/src
cmsenv
cd -
runTuple $1 $PWD/output
EOF1
chmod u+x jobs/runTuple_$FILENAME.sh
chmod u+x jobs/runTuple_$FILENAME.submit
condor_submit jobs/runTuple_$FILENAME.submit
