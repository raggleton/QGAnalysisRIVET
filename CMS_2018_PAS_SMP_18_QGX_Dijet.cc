// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Jet.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "Math/SVector.h"

#include <algorithm>


using std::cout;
using std::endl;
using std::vector;
using namespace fastjet;


namespace Rivet {

  /// @brief Class to calculate lambda variables from jet consituents
  /// Note that we ask for the daughters and 4-vector separately, to allow for
  /// custom filtering of daughters first
  class LambdaCalculator {
  public:
    LambdaCalculator(const vector<PseudoJet> & daughters, float jet_radius, const PseudoJet & jet_vector, const PseudoJet & wta_jet_vector):
      daughters_(daughters),
      jetRadius_(jet_radius),
      ptSum_(0),
      jetVector_(jet_vector),
      wtaJetVector_(wta_jet_vector)
    {
      // cache pt_sum
      for (auto & dtr : daughters_) {
        ptSum_ += dtr.pt();
      }
    }

    float getLambda(float kappa, float beta) {
      float result = 0.;
      for (auto dtr : daughters_) {
        float z = (kappa != 0) ? dtr.pt() / ptSum_ : 1.;
        // Use a dfiferent jet vector depending on beta
        float theta = (beta != 0) ? dtr.delta_R((beta <= 1) ? wtaJetVector_ : jetVector_) / jetRadius_ : 1.;
        result += (pow(z, kappa) * pow(theta, beta));
      }
      return result;
    }

  private:
    vector<PseudoJet> daughters_;
    float jetRadius_, ptSum_;
    PseudoJet jetVector_, wtaJetVector_;
  };

/**
 * Lightweight class to hold info about Lambda variable
 */
  class LambdaVar {

  public:
    LambdaVar(const std::string & name_, float kappa_, float beta_, bool isCharged_):
      name(name_),
      kappa(kappa_),
      beta(beta_),
      isCharged(isCharged_)
    {}

    std::string name;
    float kappa;
    float beta;
    bool isCharged;
  };


  /// @brief Routine for QG substructure analysis
  class CMS_2018_PAS_SMP_18_QGX_DIJET : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2018_PAS_SMP_18_QGX_DIJET);

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Particles for the jets
      FinalState fs(-5, 5, 0.0*GeV);
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      // jet_input.addVetoPairId(PID::MUON);
      addProjection(jet_input, "JET_INPUT");

      // Book histograms
      // resize vectors appropriately
      uint nHistsRadii = _jetRadii.size();
      uint nHistsLambda = _lambdaVars.size();
      uint nHistsPt = _ptBinsGen.size()-1;
      _h_dijet_cen.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
      _h_dijet_cen_groomed.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
      _h_dijet_fwd.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));
      _h_dijet_fwd_groomed.resize(nHistsRadii, vector<vector<Histo1DPtr> >(nHistsLambda, vector<Histo1DPtr>(nHistsPt)));

      // Now book histos
      // remember 1-indexed
      for (uint radiusInd=1; radiusInd <= _jetRadii.size(); radiusInd++) {
        for (uint lambdaInd=1; lambdaInd <= _lambdaVars.size(); lambdaInd++) {
          for (uint ptInd=1; ptInd < _ptBinsGen.size(); ptInd++) {
            _h_dijet_cen[radiusInd-1][lambdaInd-1][ptInd-1] = bookHisto1D((10*radiusInd) + 1, lambdaInd, ptInd);
            _h_dijet_cen_groomed[radiusInd-1][lambdaInd-1][ptInd-1] = bookHisto1D((10*radiusInd) + 2, lambdaInd, ptInd);
            _h_dijet_fwd[radiusInd-1][lambdaInd-1][ptInd-1] = bookHisto1D((10*radiusInd) + 3, lambdaInd, ptInd);
            _h_dijet_fwd_groomed[radiusInd-1][lambdaInd-1][ptInd-1] = bookHisto1D((10*radiusInd) + 4, lambdaInd, ptInd);
          }
        }
      }

      // _h_pt_jet = bookHisto1D(1, 90, 1);

      vector<double> binEdges = {
        1E-9,
        5E-9,
        1E-8,
        5E-8,
        1E-7,
        5E-7,
        1E-6,
        5E-6,
        1E-5,
        5E-5,
        1E-4,
        5E-4,
        1E-3,
        5E-3,
        1E-2,
        5E-2,
        1E-1,
        5E-1,
        1,
        5,
        10,
      };

      _h_weight = bookHisto1D("weight_jet_pt", binEdges);
      _h_jet_pt_ov_scale = bookHisto1D("jet_pt_ov_scale", 100, 0, 5);

      wta_cluster_ = Recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme),
                               false, Recluster::keep_only_hardest);
      ca_cluster_ = Recluster(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme),
                              false, Recluster::keep_only_hardest);
      mmdt_.reset(new fastjet::contrib::ModifiedMassDropTagger(0.1));
      mmdt_->set_grooming_mode();
      mmdt_->set_reclustering(false);
    }


    // Get index of largest bin smaller than value in vector
    // e.g. what you'd need when binning a continuous variable
    uint getBinIndex(float value, const vector<float> & bins) {
      auto itr = std::lower_bound(bins.begin(), bins.end(), value);
      return itr - bins.begin() - 1;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Convert Particles into PseudoJets for clustering
      const VetoedFinalState & fs = applyProjection<VetoedFinalState>(event, "JET_INPUT");
      const ParticleVector & fsParticles = fs.particles();
      vector<PseudoJet> particles;
      particles.reserve(fsParticles.size());
      for (uint iFS=0; iFS<fsParticles.size(); iFS++){
        PseudoJet p = fsParticles[iFS].pseudojet();
        p.set_user_index(iFS); // for later reference to original Particles
        particles.push_back(p);
      }

      for (uint radiusInd=0; radiusInd < _jetRadii.size(); radiusInd++) {
        float jetRadius = _jetRadii.at(radiusInd);

        JetDefinition jet_def(antikt_algorithm, jetRadius);
        double etaCut = 2.5-(jetRadius/2.);
        vector<PseudoJet> jets = (SelectorNHardest(2) * SelectorAbsRapMax(etaCut) * SelectorPtMin(15))(jet_def(particles));

        bool passDijet = false;
        if (jets.size() < 2) continue;
        const auto & jet1 = jets.at(0);
        const auto & jet2 = jets.at(1);
        float jet1pt = jet1.pt();
        float jet2pt = jet2.pt();
        float asym = (jet1pt - jet2pt) / (jet1pt+jet2pt);
        float dphi = Rivet::deltaPhi(jet1.phi(), jet2.phi());
        auto * genEvt = event.genEvent();
        // cut on this to avoid large-weight events from high jet pt, low ptHat
        float pt_ov_scale = jet1pt / genEvt->event_scale();
        passDijet = ((asym < 0.3) && (dphi > 2.0) && (pt_ov_scale < 2.));

        if (!passDijet) continue;

        // Sort by increasing absolute eta
        vector<PseudoJet> dijets = {jet1, jet2};
        std::sort(dijets.begin(), dijets.end(),
                  [] (const PseudoJet & A, const PseudoJet & B)
                  { return fabs(A.eta()) < fabs(B.eta()); }
                 );

        for (uint iJ=0; iJ<dijets.size(); iJ++) {
          bool isCentral = (iJ == 0);
          PseudoJet & jetItr = dijets[iJ];

          // _h_pt_jet->fill(jetItr.pt(), weight);

          // Simplify life - ignore this jet if it is below 1st hist pt range
          // Note that we don't apply it to the original jet pt cut - since
          // we have phase space where one jet is > 50, and one < 50
          if (jetItr.pt() < _ptBinsGen[0]) continue;
          // ignore jet if beyond the last bin
          if (jetItr.pt() > _ptBinsGen.back()) continue;

          _h_weight->fill(weight);
          _h_jet_pt_ov_scale->fill(jetItr.pt() / genEvt->event_scale());

          // Need to use original, ungroomed jet pT to bin
          uint ptBinInd = getBinIndex(jetItr.pt(), _ptBinsGen);

          // UNGROOMED VERSION
          // -------------------------------------------------------------------

          // Apply pT cut to constituents first
          vector<PseudoJet> constits = SelectorPtMin(_constitPtMin)(jetItr.constituents());

          // since this is the superset, if this fails, nothing else will pass
          if (!(constits.size() >= _minNumConstits)) continue;

          // Setup calculators
          // Get WTA jet axis used for Lambda calculations
          PseudoJet wtaJet = wta_cluster_(jetItr);
          LambdaCalculator lambdaCalc(constits, jetRadius, jetItr, wtaJet);

          // Do a version with charged-only constituents
          vector<PseudoJet> chargedConstits;
          foreach (const PseudoJet & p, constits) {
            if (fsParticles.at(p.user_index()).isCharged()) chargedConstits.push_back(p);
          }

          bool passCharged = (chargedConstits.size() >= _minNumConstits);
          LambdaCalculator lambdaCalcCharged(chargedConstits, jetRadius, jetItr, wtaJet);

          // Fill hists for each lambda variable
          for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
            const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
            if (thisLambdaVar.isCharged && !passCharged) continue;

            LambdaCalculator & thisLambdaCalc = (thisLambdaVar.isCharged) ? lambdaCalcCharged : lambdaCalc;
            float val = thisLambdaCalc.getLambda(thisLambdaVar.kappa, thisLambdaVar.beta);

            if (isCentral) {
              _h_dijet_cen[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            } else {
              _h_dijet_fwd[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            }
          }

          // GROOMED VERSION
          // -------------------------------------------------------------------
          // Get groomed jet
          PseudoJet caJet = ca_cluster_(jetItr);
          PseudoJet groomedJet = (*mmdt_)(caJet);

          // Apply pT cut to constituents first
          vector<PseudoJet> groomedConstits = SelectorPtMin(_constitPtMin)(groomedJet.constituents());
          if (!(groomedConstits.size() >= _minNumConstits)) continue;

          // Setup calculators
          // Use groomed axis instead to calculate theta
          // Get WTA jet axis used for Lambda calculations
          PseudoJet wtaGroomedJet = wta_cluster_(groomedJet);
          LambdaCalculator lambdaCalcGroomed(groomedConstits, jetRadius, groomedJet, wtaGroomedJet);

          // Do a version with charged-only constituents
          vector<PseudoJet> chargedConstitsGroomed;
          foreach (const PseudoJet & p, groomedConstits) {
            if (fsParticles.at(p.user_index()).isCharged()) chargedConstitsGroomed.push_back(p);
          }

          bool passChargedGroomed = (chargedConstitsGroomed.size() >= _minNumConstits);
          LambdaCalculator lambdaCalcChargedGroomed(chargedConstitsGroomed, jetRadius, groomedJet, wtaGroomedJet);

          // Fill hists for each lambda variable
          for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
            const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
            if (thisLambdaVar.isCharged && !passChargedGroomed) continue;

            LambdaCalculator & thisLambdaCalc = (thisLambdaVar.isCharged) ? lambdaCalcChargedGroomed : lambdaCalcGroomed;
            float val = thisLambdaCalc.getLambda(thisLambdaVar.kappa, thisLambdaVar.beta);

            if (isCentral) {
              _h_dijet_cen_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            } else {
              _h_dijet_fwd_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            }
          }

        } // end loop over dijets
      } // end loop over jet radii
    } // end analyze() function

    /**
     * Convert TH2 to TMatrixD. FIXME: errors not included. Use std::pair?
     * @param  hist    TH2* to be converted
     * @param  oflow_x Include under/overflow bins on x axis
     * @param  oflow_y Include under/overflow bins on x axis
     * @return         TMatrixD representation of TH2's bin contents.
     */
    TMatrixD th2_to_tmatrix(TH2 * hist, bool oflow_x, bool oflow_y) {
      int ncol = hist->GetNbinsX();
      if (oflow_x) ncol += 2;
      int nrow = hist->GetNbinsY();
      if (oflow_y) nrow += 2;

      TMatrixD result(nrow, ncol);

      // Get ROOT indices to loop over
      int y_start = (oflow_y) ? 0 : 1;
      int y_end = hist->GetNbinsY();
      if (oflow_y) y_end += 1;

      int x_start = (oflow_x) ? 0 : 1;
      int x_end = hist->GetNbinsX();
      if (oflow_x) x_end += 1;

      // y_ind, x_ind for indexing as always starts at 0
      // iy, ix for TH2 bin counting
      int y_ind(0), iy(y_start);
      for (; iy <= y_end; y_ind++, iy++) {
        int x_ind(0), ix(x_start);
        for (; ix <= x_end; x_ind++, ix++) {
          result(y_ind, x_ind) = hist->GetBinContent(ix, iy);
        }
      }
      return result;
    }

    /**
     * Convert YODA Histo1DPtr to TVectorD
     * @param  hist    Histo1DPtr to be converted
     * @param  oflow_x Include under/overflow bins
     * @return         TVectorD representation of Histo1D's bin contents
     */
    TVectorD histo1d_to_tvector(Histo1DPtr hist, bool oflow_x) {
      int nbins = hist->numBins();
      if (oflow_x) nbins += 2;
      cout << "Creating vector with " << nbins << " bins" << endl;
      TVectorD result(nbins);

      // Fill under/overflow into vector
      if (oflow_x) {
        cout << "Setting bin 0" << endl;
        result[0] = hist->underflow().sumW();
        cout << "Setting bin " << nbins-1 << endl;
        result[nbins-1] = hist->overflow().sumW();
      }

      // x_ind for vector index
      // ix for Yoda bin counting
      int x_ind = (oflow_x) ? 1 : 0;
      uint ix(0);
      cout << "histo1d_to_tvector:" << endl;
      for (; ix < hist->numBins(); ix++, x_ind++) {
        cout << "Yoda Bin " << ix << " sumW: " << hist->bin(ix).sumW() << " area: " << hist->bin(ix).area() << endl;
        result[x_ind] = hist->bin(ix).sumW();
      }

      return result;
    }

    /**
     * Normalise TMatrix such that sum over a given column = 1
     */
    void normalise_tmatrix_by_col(TMatrixD & matrix) {
      int nrow = matrix.GetNrows();
      int ncol = matrix.GetNcols();
      cout << "normalise_tmatrix_by_col" << endl;
      cout << "nrow: " << nrow << " ncol: " << ncol << endl;

      // create vector of weights, one entry per column
      TVectorD weights(ncol);
      for (int i=0; i < ncol; i++) {
        cout << "Doing col " << i << endl;
        // Get bin contents for this col
        TMatrixDColumn thisCol(matrix, i);

        float this_col_sum = 0;
        for (int j=0; j < nrow; j++) {
          this_col_sum += thisCol(j);
        }
        cout << "this_col_sum " << this_col_sum << endl;

        if (this_col_sum != 0) {
          weights[i] = 1./this_col_sum;
        } else {
          weights[i] = 1.;
        }
      }

      cout << "Got weights, now normalise" << endl;
      matrix.NormByRow(weights, "");  // MUST use empty string otherwise weights applied as inverse
      cout << "Done normalising" << endl;
    }

    /**
     * Do folding of vector by response_matrix. Automatically normalises matrix
     * by column (gen bins)
     */
    TVectorD get_folded_tvector(TMatrixD & response_matrix, TVectorD & vector) {

      // Normalise response_matrix so that bins represent prob to go from
      // given gen bin to a reco bin
      // ASSUMES GEN ON X AXIS!
      normalise_tmatrix_by_col(response_matrix);

      // Multiply
      cout << "Multiplying" << endl;
      cout << "vector size: " << vector.GetNrows() << endl;
      TVectorD folded_vec = response_matrix * vector;
      cout << "done Multiplying" << endl;

      return folded_vec;
    }

    /**
     * Place contents for Histo1DPtr with that in vector
     * @param oflow_x If true, then 1st and last entries in vector correspond to under/overflow bins
     */
    void refil_histo1d_from_tvector(Histo1DPtr hist, TVectorD & vector, bool oflow_x) {

      if (oflow_x) {
        // Do manual filling for under/overflow
        hist->underflow().reset();
        hist->underflow().fill(vector[0]);
        hist->overflow().reset();
        hist->overflow().fill(vector[vector.GetNrows()-1]);
      }

      // x_ind for vector index
      // ix for Yoda bin counting
      int x_ind = (oflow_x) ? 1 : 0;
      uint ix(0);
      // cout << "refil_histo1d_from_tvector:" << endl;
      for (; ix < hist->numBins(); ix++) {
        // cout << "Bin " << ix << " = " << vector[x_ind] << endl;
        hist->fill(hist->bin(ix).xMid(), vector[x_ind]);
        cout << hist->bin(ix).area() << endl;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      float lumi = 35918;

      // norm hists to cross section
      // scale(_h_pt_jet, lumi * crossSection() / sumOfWeights());

      // normalise hists to unity
      // TODO does it divide by bin width?
      // for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
      //   for (uint ptInd=0; ptInd < _ptBinsGen.size()-1; ptInd++) {
      //     normalize(_h_dijet_cen[lambdaInd][ptInd]);
      //   }
      // }

      // TFile * _response_file = TFile::Open("uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root");
      // // TFile * _response_file = TFile::Open("pythia_response.root");

      // // Get response matrix
      // TH2D * _response_hist = (TH2D*) _response_file->Get("Dijet_QG_Unfold_tighter/tu_pt_GenReco_all");
      // // rebin the detector axis as half-binned
      // _response_hist->RebinY(2);

      // // Do folding & updating of Histo1D
      // bool oflow = true;
      // TMatrixD matrix_response = th2_to_tmatrix(_response_hist, oflow, oflow);
      // cout << "MATRIX SIZE: " << matrix_response.GetNrows() << " x " << matrix_response.GetNcols() << endl;
      // // cout << "Matrix(5,5) = " << matrix_response(5, 5) << endl;
      // TVectorD pt_vector = histo1d_to_tvector(_h_pt_jet, oflow);
      // cout << "VECTOR SIZE: " << pt_vector.GetNrows() << endl;
      // TVectorD pt_folded = get_folded_tvector(matrix_response, pt_vector);
      // // cout << "Matrix(5,5) normed = " << matrix_response(5, 5) << endl;
      // _h_pt_jet->reset();
      // refil_histo1d_from_tvector(_h_pt_jet, pt_folded, oflow);

    } // end of finalize

    // Order matters here
    const vector<float> _jetRadii = {0.4, 0.8};

    // This order is important! index in vector used to create YODA plot name
    // Must match that in extracRivetPlotsDijet.py
    const vector<LambdaVar> _lambdaVars = {
      LambdaVar("jet_puppiMultiplicity", 0, 0, false),
      LambdaVar("jet_pTD", 2, 0, false),
      LambdaVar("jet_LHA", 1, 0.5, false),
      LambdaVar("jet_width", 1, 1, false),
      LambdaVar("jet_thrust", 1, 2, false),
      LambdaVar("jet_puppiMultiplicity_charged", 0, 0, true),
      LambdaVar("jet_pTD_charged", 2, 0, true),
      LambdaVar("jet_LHA_charged", 1, 0.5, true),
      LambdaVar("jet_width_charged", 1, 1, true),
      LambdaVar("jet_thrust_charged", 1, 2, true),
    };

    const vector<float> _ptBinsGen = {
      50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000
    };

    /// @name Histograms
    Histo1DPtr _h_pt_jet;
    Histo1DPtr _h_weight, _h_jet_pt_ov_scale;

    // 3D vector: [jet radius][lambda variable][pt bin]
    // since each pt bin has its own normalised distribution
    vector<vector<vector<Histo1DPtr> > > _h_dijet_cen,
                                         _h_dijet_cen_groomed,
                                         _h_dijet_fwd,
                                         _h_dijet_fwd_groomed;

    const float _constitPtMin = 1.;
    const uint _minNumConstits = 2;

    Recluster wta_cluster_;
    Recluster ca_cluster_;
    std::unique_ptr<fastjet::contrib::ModifiedMassDropTagger> mmdt_;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_PAS_SMP_18_QGX_DIJET);


}
