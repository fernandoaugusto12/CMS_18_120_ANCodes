#!/bin/bash

melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc530
melaincdir=${CMSSW_BASE}/src


cmsswlibdir=$CMSSW_RELEASE_BASE/lib/slc6_amd64_gcc530
cmsswincdir=$CMSSW_RELEASE_BASE/src

export LD_LIBRARY_PATH=${melalibdir}:${cmsswlibdir}:$LD_LIBRARY_PATH


g++ -I $ROOTSYS/include -I $ROOFITSYS/include -I ${melaincdir} -I ${cmsswincdir} ComputeFRandCRs.cc `root-config --glibs` `root-config --libs` `root-config --cflags` -L $ROOFITSYS/lib  -lRooFit -lRooFitCore -L ${melalibdir} -L ${cmsswlibdir}  -lZZMatrixElementMELA   -lZZMatrixElementMEMCalculators -lKaMuCaCalibration -lCondFormatsJetMETObjects -l JetMETCorrectionsModules  -o ComputeFRandCRs
