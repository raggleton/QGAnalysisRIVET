# CMS Quark-Gluon analysis Rivet routines

There are 2 routines: one for dijet selection ([`CMS_2018_PAS_SMP_18_QGX_Dijet.cc`](CMS_2018_PAS_SMP_18_QGX_Dijet.cc)), and one for Z+jet selection ([`CMS_2018_PAS_SMP_18_QGX_ZPJ.cc`](CMS_2018_PAS_SMP_18_QGX_ZPJ.cc)).

There are also various corresponding `.yoda`, `.plot`, and CMSSW python configs, explained below.
The yoda files contain all the histograms at particle level, and were produced in the main analysis using a particular generator configuration.
Each routine has a corresponding `.plot` file to correctly label the histograms.


## Dijet routine

This is designed to be used with QCD multijet event generation.

The phase space is divided into:

- separate central & forward jet histograms (from the dijet selection)
- AK4 & AK8 jets
- charged+neutral, and charged-only lambda variables (the 5 specified in 1704.03878)
- ungroomed & groomed versions of the variables (groomed with Modified Mass-Drop Tagging)
- separate jet pT bins

An example CMSSW routine with the generator configuration is in [`RIVET_QCD_Pt-15To7000_TuneCUETP8M1_Flat_13TeV-pythia8_cff.py`](RIVET_QCD_Pt-15To7000_TuneCUETP8M1_Flat_13TeV-pythia8_cff.py), generating PYTHIA8 QCD events with the P8M1 tune.

A few example YODA files have been included, produced by different generators:

- [`CMS_2018_PAS_SMP_18_QGX_DIJET_GEN_MGPYTHIA.yoda`](CMS_2018_PAS_SMP_18_QGX_DIJET_GEN_MGPYTHIA.yoda): MadGraph5_2.2.2 + PYTHIA 8.212, P8M1 tune
- [`CMS_2018_PAS_SMP_18_QGX_DIJET_GEN_PYTHIA_FLAT.yoda`](CMS_2018_PAS_SMP_18_QGX_DIJET_GEN_PYTHIA_FLAT.yoda): PYTHIA 8.?, P8M1 tune


### Histogram naming scheme

YODA Histograms have names `d<DD>-x<XX>-y<YY>`:

- `<DD>`: "channel" =  {AK4, AK8} [10, 20] + {dijet central, dijet central groomed, dijet forward, dijet forward groomed} [01, ..., 04]. e.g. 22 is central AK8 groomed jet.
- `<XX>`: lambda variable, 01-05 are neutral+charged, 06-10 are charged-only
- `<YY>`: jet pT bin


## Z+jet routine

This is designed to be used with Z+jet event generation.

The phase space is divided into:

- AK4 & AK8 jet (note we only consider the 1st jet)
- charged+neutral, and charged-only lambda variables (the 5 specified in 1704.03878)
- ungroomed & groomed versions of the variables (groomed with Modified Mass-Drop Tagging)
- separate jet pT bins

An example YODA file has been included, [`CMS_2018_PAS_SMP_18_QGX_ZPJ_GEN_MGPYTHIA.yoda`](CMS_2018_PAS_SMP_18_QGX_ZPJ_GEN_MGPYTHIA.yoda), made by MG5 + Pythia8 with the P8M1 tune.

### Histogram naming scheme

YODA Histograms have names `d<DD>-x<XX>-y<YY>`:

- `<DD>`: "channel" = {AK4, AK8} [10, 20] + {ungroomed, groomed} [01, 02]. e.g. 21 is AK8 jet, ungroomed.
- `<XX>`: lambda variable, 01-05 are neutral+charged, 06-10 are charged-only
- `<YY>`: jet pT bin

