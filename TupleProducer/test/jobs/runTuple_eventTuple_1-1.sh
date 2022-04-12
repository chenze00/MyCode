source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
#cd /afs/cern.ch/user/c/chenz/CMSSW_10_6_20/src
cd /afs/cern.ch/user/c/chenz/CMSSW_10_6_20/src
cmsenv
cd -
runTuple /eos/cms/store/group/phys_tau/TauML/prod_2018_v2/full_tuples/DYJetsToLL_M-50/eventTuple_1-1.root output
