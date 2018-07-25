#!/bin/bash -e

PRE=../CMSSW_10_0_0/src

rsync $PRE/Rivet/SMP/data/CMS_2018_PAS_SMP_18_QGX.plot .
rsync $PRE/Rivet/SMP/data/CMS_2018_PAS_SMP_18_QGX.info .
rsync $PRE/Rivet/SMP/src/CMS_2018_PAS_SMP_18_QGX.cc .

rsync $PRE/root2yoda .
rsync $PRE/RIVET_QCD_Pt-15To7000_TuneCUETP8M1_Flat_13TeV-pythia8_cff.py .


