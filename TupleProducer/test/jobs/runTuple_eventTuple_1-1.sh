source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cd /eos/home-c/chenz/MyCode/CMSSW_10_6_20/src
cmsenv
cd -
runTuple eventTuple_1-1 output
