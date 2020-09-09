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
#include "fastjet/contrib/SoftDrop.hh"
//#include "fastjet/contrib/ModifiedMassDropTagger.hh"

#include <algorithm>


using std::cout;
using std::endl;
using std::vector;
using namespace fastjet;


namespace Rivet {

/// \class Angularity
/// definition of angularity
///
class Angularity : public FunctionOfPseudoJet<double>{
public:
  /// ctor
  Angularity(double alpha, double jet_radius, double kappa=1.0, bool isCharged=false, Selector constitCut=SelectorPtMin(0.)) : _alpha(alpha), _radius(jet_radius), _kappa(kappa), _isCharged(isCharged), _constitCut(constitCut) {}

  /// description
  std::string description() const{
    ostringstream oss;
    oss << "Angularity with alpha=" << _alpha;
    return oss.str();
  }

  /// computation of the angularity itself
  double result(const PseudoJet &jet) const{
    // get the jet constituents
    vector<PseudoJet> constits = jet.constituents();

    // get the reference axis
    PseudoJet reference_axis = _get_reference_axis(jet);

    // do the actual coputation
    double numerator = 0.0, denominator = 0.0;
    unsigned int num = 0;
    for (const auto &c : constits){
      if ((_isCharged) && (!c.user_index())) continue;
      if (!_constitCut.pass(c)) continue;
      double pt = c.pt();
      // Note: better compute (dist^2)^(alpha/2) to avoid an extra square root
      numerator   += pow(pt, _kappa) * pow(c.squared_distance(reference_axis), 0.5*_alpha);
      denominator += pt;
      num += 1;
    }
    if (!((num >= _minNumConstits) && (denominator > 0))) return -1;
    // the formula is only correct for the the typical angularities which satisfy either kappa==1 or beta==0.
    else return numerator/(pow(denominator, _kappa)*pow(_radius, _alpha));
  }

protected:
  PseudoJet _get_reference_axis(const PseudoJet &jet) const{
    if (_alpha>1) return jet;

    Recluster recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme));
    return recluster(jet);
  }

  double _alpha, _radius, _kappa;
  bool _isCharged;
  Selector _constitCut;
  const uint _minNumConstits = 2;
};

/**
 * Lightweight class to hold info about Lambda variable
 */
  class LambdaVar {

  public:
    LambdaVar(const std::string & name_, float kappa_, float beta_, bool isCharged_, Selector constitCut_):
      name(name_),
      kappa(kappa_),
      beta(beta_),
      isCharged(isCharged_),
      constitCut(constitCut_)
    {}

    std::string name;
    float kappa;
    float beta;
    bool isCharged;
    Selector constitCut;
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

      // yoda plot naming scheme
      // --------------------------------------------------------------------------
      // d = channel (ak4/8 radii (10, 20) + {dijet cen / fwd} * groomed versions[1..4])
      // x = lambda variable; neutral+charged & charged-only are treated separately
      // y = pT bin
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
        p.set_user_index(fsParticles[iFS].isCharged()); // for later reference to charge
        particles.push_back(p);
      }

      for (uint radiusInd=0; radiusInd < _jetRadii.size(); radiusInd++) {
        float jetRadius = _jetRadii.at(radiusInd);

        JetDefinition jet_def(antikt_algorithm, jetRadius);
        vector<PseudoJet> jets = (SelectorNHardest(2) * SelectorAbsRapMax(1.7) * SelectorPtMin(15))(jet_def(particles));

        bool passDijet = false;
        if (jets.size() < 2) continue;
        const auto & jet1 = jets.at(0);
        const auto & jet2 = jets.at(1);
        float jet1pt = jet1.pt();
        float jet2pt = jet2.pt();
        float asym = (jet1pt - jet2pt) / (jet1pt+jet2pt);
        float dphi = Rivet::deltaPhi(jet1.phi(), jet2.phi());
        passDijet = ((asym < 0.3) && (dphi > 2.0));

        if (!passDijet) continue;

        // Sort by increasing absolute rapidity
        vector<PseudoJet> dijets = {jet1, jet2};
        std::sort(dijets.begin(), dijets.end(),
                  [] (const PseudoJet & A, const PseudoJet & B)
                  { return fabs(A.rapidity()) < fabs(B.rapidity()); }
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

          // Need to use original, ungroomed jet pT to bin
          uint ptBinInd = getBinIndex(jetItr.pt(), _ptBinsGen);

          // UNGROOMED VERSION
          // -------------------------------------------------------------------

          // Fill hists for each lambda variable
          for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
            const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
            Angularity angularity(thisLambdaVar.beta, jetRadius, thisLambdaVar.kappa, thisLambdaVar.isCharged, thisLambdaVar.constitCut);
            float val = angularity(jetItr);
            if (val<0) continue;

            if (isCentral) {
              _h_dijet_cen[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            } else {
              _h_dijet_fwd[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            }
          }

          // GROOMED VERSION
          // -------------------------------------------------------------------
          // Get groomed jet
          fastjet::contrib::SoftDrop sd(0, 0.1, jetRadius);
          PseudoJet groomedJet = sd(jetItr);
          //fastjet::contrib::ModifiedMassDropTagger mmdt(0.1); mmdt.set_grooming_mode(); mmdt.set_reclustering(false);
          //Recluster ca_cluster(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R));
          //PseudoJet groomedJet = mmdt(ca_cluster(jetItr));


          // Fill hists for each lambda variable
          for (uint lambdaInd=0; lambdaInd < _lambdaVars.size(); lambdaInd++) {
            const LambdaVar & thisLambdaVar = _lambdaVars[lambdaInd];
            Angularity angularity(thisLambdaVar.beta, jetRadius, thisLambdaVar.kappa, thisLambdaVar.isCharged, thisLambdaVar.constitCut);
            float val = angularity(groomedJet);
            if (val<0) continue;

            if (isCentral) {
              _h_dijet_cen_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            } else {
              _h_dijet_fwd_groomed[radiusInd][lambdaInd][ptBinInd]->fill(val, weight);
            }
          }

        } // end loop over dijets
      } // end loop over jet radii
    } // end analyze() function


    /// Normalise histograms etc., after the run
    void finalize() {

      // norm hists to cross section
      // scale(_h_pt_jet, lumi * crossSection() / sumOfWeights());

    } // end of finalize

    // Order matters here
    const vector<float> _jetRadii = {0.4, 0.8};

    // This order is important! index in vector used to create YODA plot name
    // Must match that in extracRivetPlotsDijet.py
    const vector<LambdaVar> _lambdaVars = {
      LambdaVar("jet_multiplicity", 0, 0, false, SelectorPtMin(1.)),
      LambdaVar("jet_pTD", 2, 0, false, SelectorPtMin(0.)),
      LambdaVar("jet_LHA", 1, 0.5, false, SelectorPtMin(0.)),
      LambdaVar("jet_width", 1, 1, false, SelectorPtMin(0.)),
      LambdaVar("jet_thrust", 1, 2, false, SelectorPtMin(0.)),
      LambdaVar("jet_multiplicity_charged", 0, 0, true, SelectorPtMin(1.)),
      LambdaVar("jet_pTD_charged", 2, 0, true, SelectorPtMin(0.)),
      LambdaVar("jet_LHA_charged", 1, 0.5, true, SelectorPtMin(0.)),
      LambdaVar("jet_width_charged", 1, 1, true, SelectorPtMin(0.)),
      LambdaVar("jet_thrust_charged", 1, 2, true, SelectorPtMin(0.)),
    };

    const vector<float> _ptBinsGen = {
      50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000
    };

    /// @name Histograms
    Histo1DPtr _h_pt_jet;

    // 3D vector: [jet radius][lambda variable][pt bin]
    // since each pt bin has its own normalised distribution
    vector<vector<vector<Histo1DPtr> > > _h_dijet_cen,
                                         _h_dijet_cen_groomed,
                                         _h_dijet_fwd,
                                         _h_dijet_fwd_groomed;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_PAS_SMP_18_QGX_DIJET);


}
