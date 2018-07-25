// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Routine for QG substructure analysis
  class CMS_2018_PAS_SMP_18_QGX : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2018_PAS_SMP_18_QGX);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs(Cuts::abseta < 5);
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      _h_pt_jet = bookHisto1D(1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);
      
      if (jets.size() >= 2) {
        float avePt = 0.5*(jets.at(0).pT() + jets.at(1).pT());
        _h_pt_jet->fill(avePt, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_pt_jet, crossSection()/picobarn/sumOfWeights()); // norm to cross section

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pt_jet;

    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_PAS_SMP_18_QGX);


}
