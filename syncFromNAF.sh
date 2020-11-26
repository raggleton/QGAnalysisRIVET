#!/bin/bash

filenames=(SMP/data/CMS_2018_PAS_SMP_18_QGX_DIJET_GEN_MGPYTHIA.yoda SMP/data/CMS_2018_PAS_SMP_18_QGX_ZPJ_GEN_MGPYTHIA.yoda SMP/data/CMS_2018_PAS_SMP_18_QGX_DIJET.plot SMP/data/CMS_2018_PAS_SMP_18_QGX_ZPJ.plot SMP/src/CMS_2018_PAS_SMP_18_QGX_ZPJ.cc SMP/src/CMS_2018_PAS_SMP_18_QGX_DIJET.cc)

nafdir="/nfs/dust/cms/user/aggleton/QG/Rivet_10_6/CMSSW_10_6_0/src/Rivet/"
for fname in ${filenames[@]}; do
    echo $fname
    rsync -avzhP NAFEL7:"$nafdir/$fname" .
done

