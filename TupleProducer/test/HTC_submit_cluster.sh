#!/bin/sh 
# $1 - full path input file

#/eos/cms/store/group/phys_tau/TauML/prod_2018_v2/full_tuples/DYJetsToLL_M-50

cat > jobs/runTuple_cluster.submit <<EOF
#!/bin/sh
+RequestRuntime = 1200

executable = jobs/runTuple_cluster.sh
arguments  = \$(file)

transfer_executable = True
universe            = vanilla
getenv              = True
Requirements        = OpSysAndVer == "CentOS7"

output              = jobs/\$(name).out
error               = jobs/\$(name).error
log                 = jobs/\$(name).log

queue 1 file,name from jobs/items.txt

EOF

cat > jobs/runTuple_cluster.sh <<EOF1
#!/bin/sh
export PYTHONHOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/python/2.7.14-pafccj 
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700  
cd ${CMSSW_BASE}/src
cmsenv
cd -
runTuple \$1 $PWD/output
EOF1
chmod u+x jobs/runTuple_cluster.sh
chmod u+x jobs/runTuple_cluster.submit
condor_submit jobs/runTuple_cluster.submit
