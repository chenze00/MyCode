//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Aug 22 15:24:28 2021 by ROOT version 6.22/08
// from TTree taus/taus
// found on file: /nfs/dust/cms/user/cardinia/TauReco/example_rootfiles/DYJetsToLL_M-50.root
//////////////////////////////////////////////////////////

#ifndef MyTauClass_h
#define MyTauClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
//#include "vector"
using namespace std;

class MyTauClass {
public :
  TString         inputName; //file name
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       evt;
   Int_t           npv;
   Float_t         rho;
   Float_t         genEventWeight;
   Float_t         trainingWeight;
   Int_t           sampleType;
   Int_t           tauType;
   ULong64_t       dataset_id;
   ULong64_t       dataset_group_id;
   Float_t         npu;
   Float_t         pv_x;
   Float_t         pv_y;
   Float_t         pv_z;
   Float_t         pv_t;
   Float_t         pv_xE;
   Float_t         pv_yE;
   Float_t         pv_zE;
   Float_t         pv_tE;
   Float_t         pv_chi2;
   Float_t         pv_ndof;
   Int_t           entry_index;
   Int_t           total_entries;
   Int_t           genLepton_index;
   Int_t           genLepton_kind;
   Int_t           genLepton_charge;
   Float_t         genLepton_vis_pt;
   Float_t         genLepton_vis_eta;
   Float_t         genLepton_vis_phi;
   Float_t         genLepton_vis_mass;
   Int_t           genLepton_lastMotherIndex;
   vector<int>     *genParticle_pdgId;
   vector<Long64_t> *genParticle_mother;
   vector<int>     *genParticle_charge;
   vector<int>     *genParticle_isFirstCopy;
   vector<int>     *genParticle_isLastCopy;
   vector<float>   *genParticle_pt;
   vector<float>   *genParticle_eta;
   vector<float>   *genParticle_phi;
   vector<float>   *genParticle_mass;
   vector<float>   *genParticle_vtx_x;
   vector<float>   *genParticle_vtx_y;
   vector<float>   *genParticle_vtx_z;
   Int_t           genJet_index;
   Float_t         genJet_pt;
   Float_t         genJet_eta;
   Float_t         genJet_phi;
   Float_t         genJet_mass;
   Float_t         genJet_emEnergy;
   Float_t         genJet_hadEnergy;
   Float_t         genJet_invisibleEnergy;
   Float_t         genJet_auxiliaryEnergy;
   Float_t         genJet_chargedHadronEnergy;
   Float_t         genJet_neutralHadronEnergy;
   Float_t         genJet_chargedEmEnergy;
   Float_t         genJet_neutralEmEnergy;
   Float_t         genJet_muonEnergy;
   Int_t           genJet_chargedHadronMultiplicity;
   Int_t           genJet_neutralHadronMultiplicity;
   Int_t           genJet_chargedEmMultiplicity;
   Int_t           genJet_neutralEmMultiplicity;
   Int_t           genJet_muonMultiplicity;
   Int_t           genJet_n_bHadrons;
   Int_t           genJet_n_cHadrons;
   Int_t           genJet_n_partons;
   Int_t           genJet_n_leptons;
   Int_t           genJet_hadronFlavour;
   Int_t           genJet_partonFlavour;
   Int_t           jet_index;
   Int_t           fatJet_index;
   Float_t         jet_pt;
   Float_t         fatJet_pt;
   Float_t         jet_eta;
   Float_t         fatJet_eta;
   Float_t         jet_phi;
   Float_t         fatJet_phi;
   Float_t         jet_mass;
   Float_t         fatJet_mass;
   Float_t         jet_neutralHadronEnergyFraction;
   Float_t         fatJet_neutralHadronEnergyFraction;
   Float_t         jet_neutralEmEnergyFraction;
   Float_t         fatJet_neutralEmEnergyFraction;
   Int_t           jet_nConstituents;
   Int_t           fatJet_nConstituents;
   Int_t           jet_chargedMultiplicity;
   Int_t           fatJet_chargedMultiplicity;
   Int_t           jet_neutralMultiplicity;
   Int_t           fatJet_neutralMultiplicity;
   Int_t           jet_partonFlavour;
   Int_t           fatJet_partonFlavour;
   Int_t           jet_hadronFlavour;
   Int_t           fatJet_hadronFlavour;
   Float_t         jet_m_softDrop;
   Float_t         fatJet_m_softDrop;
   Float_t         jet_nJettiness_tau1;
   Float_t         fatJet_nJettiness_tau1;
   Float_t         jet_nJettiness_tau2;
   Float_t         fatJet_nJettiness_tau2;
   Float_t         jet_nJettiness_tau3;
   Float_t         fatJet_nJettiness_tau3;
   Float_t         jet_nJettiness_tau4;
   Float_t         fatJet_nJettiness_tau4;
   vector<float>   *jet_subJet_pt;
   vector<float>   *fatJet_subJet_pt;
   vector<float>   *jet_subJet_eta;
   vector<float>   *fatJet_subJet_eta;
   vector<float>   *jet_subJet_phi;
   vector<float>   *fatJet_subJet_phi;
   vector<float>   *jet_subJet_mass;
   vector<float>   *fatJet_subJet_mass;
   Int_t           tau_index;
   Int_t           boostedTau_index;
   Float_t         tau_pt;
   Float_t         boostedTau_pt;
   Float_t         tau_eta;
   Float_t         boostedTau_eta;
   Float_t         tau_phi;
   Float_t         boostedTau_phi;
   Float_t         tau_mass;
   Float_t         boostedTau_mass;
   Int_t           tau_charge;
   Int_t           boostedTau_charge;
   Int_t           tau_decayMode;
   Int_t           boostedTau_decayMode;
   Int_t           tau_decayModeFinding;
   Int_t           boostedTau_decayModeFinding;
   Int_t           tau_decayModeFindingNewDMs;
   Int_t           boostedTau_decayModeFindingNewDMs;
   Float_t         tau_chargedIsoPtSum;
   Float_t         boostedTau_chargedIsoPtSum;
   Float_t         tau_chargedIsoPtSumdR03;
   Float_t         boostedTau_chargedIsoPtSumdR03;
   Float_t         tau_footprintCorrection;
   Float_t         boostedTau_footprintCorrection;
   Float_t         tau_footprintCorrectiondR03;
   Float_t         boostedTau_footprintCorrectiondR03;
   Float_t         tau_neutralIsoPtSum;
   Float_t         boostedTau_neutralIsoPtSum;
   Float_t         tau_neutralIsoPtSumWeight;
   Float_t         boostedTau_neutralIsoPtSumWeight;
   Float_t         tau_neutralIsoPtSumWeightdR03;
   Float_t         boostedTau_neutralIsoPtSumWeightdR03;
   Float_t         tau_neutralIsoPtSumdR03;
   Float_t         boostedTau_neutralIsoPtSumdR03;
   Float_t         tau_photonPtSumOutsideSignalCone;
   Float_t         boostedTau_photonPtSumOutsideSignalCone;
   Float_t         tau_photonPtSumOutsideSignalConedR03;
   Float_t         boostedTau_photonPtSumOutsideSignalConedR03;
   Float_t         tau_puCorrPtSum;
   Float_t         boostedTau_puCorrPtSum;
   UShort_t        tau_byCombinedIsolationDeltaBetaCorr3Hits;
   UShort_t        boostedTau_byCombinedIsolationDeltaBetaCorr3Hits;
   Float_t         tau_byCombinedIsolationDeltaBetaCorr3Hitsraw;
   Float_t         boostedTau_byCombinedIsolationDeltaBetaCorr3Hitsraw;
   UShort_t        tau_byDeepTau2017v2p1VSe;
   UShort_t        boostedTau_byDeepTau2017v2p1VSe;
   Float_t         tau_byDeepTau2017v2p1VSeraw;
   Float_t         boostedTau_byDeepTau2017v2p1VSeraw;
   UShort_t        tau_byDeepTau2017v2p1VSmu;
   UShort_t        boostedTau_byDeepTau2017v2p1VSmu;
   Float_t         tau_byDeepTau2017v2p1VSmuraw;
   Float_t         boostedTau_byDeepTau2017v2p1VSmuraw;
   UShort_t        tau_byDeepTau2017v2p1VSjet;
   UShort_t        boostedTau_byDeepTau2017v2p1VSjet;
   Float_t         tau_byDeepTau2017v2p1VSjetraw;
   Float_t         boostedTau_byDeepTau2017v2p1VSjetraw;
   UShort_t        tau_byIsolationMVArun2017v2DBoldDMwLT2017;
   UShort_t        boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017;
   Float_t         tau_byIsolationMVArun2017v2DBoldDMwLT2017raw;
   Float_t         boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017raw;
   UShort_t        tau_byIsolationMVArun2017v2DBnewDMwLT2017;
   UShort_t        boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017;
   Float_t         tau_byIsolationMVArun2017v2DBnewDMwLT2017raw;
   Float_t         boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017raw;
   UShort_t        tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   UShort_t        boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017;
   Float_t         tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw;
   Float_t         boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw;
   UShort_t        tau_byIsolationMVADBnewDMwLTPhase2;
   UShort_t        boostedTau_byIsolationMVADBnewDMwLTPhase2;
   Float_t         tau_byIsolationMVADBnewDMwLTPhase2raw;
   Float_t         boostedTau_byIsolationMVADBnewDMwLTPhase2raw;
   Float_t         tau_dxy_pca_x;
   Float_t         boostedTau_dxy_pca_x;
   Float_t         tau_dxy_pca_y;
   Float_t         boostedTau_dxy_pca_y;
   Float_t         tau_dxy_pca_z;
   Float_t         boostedTau_dxy_pca_z;
   Float_t         tau_dxy;
   Float_t         boostedTau_dxy;
   Float_t         tau_dxy_error;
   Float_t         boostedTau_dxy_error;
   Float_t         tau_ip3d;
   Float_t         boostedTau_ip3d;
   Float_t         tau_ip3d_error;
   Float_t         boostedTau_ip3d_error;
   Float_t         tau_dz;
   Float_t         boostedTau_dz;
   Float_t         tau_dz_error;
   Float_t         boostedTau_dz_error;
   Int_t           tau_hasSecondaryVertex;
   Int_t           boostedTau_hasSecondaryVertex;
   Float_t         tau_sv_x;
   Float_t         boostedTau_sv_x;
   Float_t         tau_sv_y;
   Float_t         boostedTau_sv_y;
   Float_t         tau_sv_z;
   Float_t         boostedTau_sv_z;
   Float_t         tau_flightLength_x;
   Float_t         boostedTau_flightLength_x;
   Float_t         tau_flightLength_y;
   Float_t         boostedTau_flightLength_y;
   Float_t         tau_flightLength_z;
   Float_t         boostedTau_flightLength_z;
   Float_t         tau_flightLength_sig;
   Float_t         boostedTau_flightLength_sig;
   Float_t         tau_pt_weighted_deta_strip;
   Float_t         boostedTau_pt_weighted_deta_strip;
   Float_t         tau_pt_weighted_dphi_strip;
   Float_t         boostedTau_pt_weighted_dphi_strip;
   Float_t         tau_pt_weighted_dr_signal;
   Float_t         boostedTau_pt_weighted_dr_signal;
   Float_t         tau_pt_weighted_dr_iso;
   Float_t         boostedTau_pt_weighted_dr_iso;
   Float_t         tau_leadingTrackNormChi2;
   Float_t         boostedTau_leadingTrackNormChi2;
   Float_t         tau_e_ratio;
   Float_t         boostedTau_e_ratio;
   Float_t         tau_gj_angle_diff;
   Float_t         boostedTau_gj_angle_diff;
   Int_t           tau_n_photons;
   Int_t           boostedTau_n_photons;
   Float_t         tau_emFraction;
   Float_t         boostedTau_emFraction;
   Int_t           tau_inside_ecal_crack;
   Int_t           boostedTau_inside_ecal_crack;
   Float_t         tau_leadChargedCand_etaAtEcalEntrance;
   Float_t         boostedTau_leadChargedCand_etaAtEcalEntrance;
   vector<int>     *pfCand_index;
   vector<int>     *lostTrack_index;
   vector<int>     *pfCand_tauSignal;
   vector<int>     *lostTrack_tauSignal;
   vector<int>     *pfCand_tauLeadChargedHadrCand;
   vector<int>     *lostTrack_tauLeadChargedHadrCand;
   vector<int>     *pfCand_tauIso;
   vector<int>     *lostTrack_tauIso;
   vector<int>     *pfCand_boostedTauSignal;
   vector<int>     *lostTrack_boostedTauSignal;
   vector<int>     *pfCand_boostedTauLeadChargedHadrCand;
   vector<int>     *lostTrack_boostedTauLeadChargedHadrCand;
   vector<int>     *pfCand_boostedTauIso;
   vector<int>     *lostTrack_boostedTauIso;
   vector<int>     *pfCand_jetDaughter;
   vector<int>     *lostTrack_jetDaughter;
   vector<int>     *pfCand_fatJetDaughter;
   vector<int>     *lostTrack_fatJetDaughter;
   vector<int>     *pfCand_subJetDaughter;
   vector<int>     *lostTrack_subJetDaughter;
   vector<float>   *pfCand_pt;
   vector<float>   *lostTrack_pt;
   vector<float>   *pfCand_eta;
   vector<float>   *lostTrack_eta;
   vector<float>   *pfCand_phi;
   vector<float>   *lostTrack_phi;
   vector<float>   *pfCand_mass;
   vector<float>   *lostTrack_mass;
   vector<int>     *pfCand_pvAssociationQuality;
   vector<int>     *lostTrack_pvAssociationQuality;
   vector<int>     *pfCand_fromPV;
   vector<int>     *lostTrack_fromPV;
   vector<float>   *pfCand_puppiWeight;
   vector<float>   *lostTrack_puppiWeight;
   vector<float>   *pfCand_puppiWeightNoLep;
   vector<float>   *lostTrack_puppiWeightNoLep;
   vector<int>     *pfCand_particleType;
   vector<int>     *lostTrack_particleType;
   vector<int>     *pfCand_charge;
   vector<int>     *lostTrack_charge;
   vector<int>     *pfCand_lostInnerHits;
   vector<int>     *lostTrack_lostInnerHits;
   vector<int>     *pfCand_nHits;
   vector<int>     *lostTrack_nHits;
   vector<int>     *pfCand_nPixelHits;
   vector<int>     *lostTrack_nPixelHits;
   vector<int>     *pfCand_nPixelLayers;
   vector<int>     *lostTrack_nPixelLayers;
   vector<int>     *pfCand_nStripLayers;
   vector<int>     *lostTrack_nStripLayers;
   vector<float>   *pfCand_vertex_x;
   vector<float>   *lostTrack_vertex_x;
   vector<float>   *pfCand_vertex_y;
   vector<float>   *lostTrack_vertex_y;
   vector<float>   *pfCand_vertex_z;
   vector<float>   *lostTrack_vertex_z;
   vector<float>   *pfCand_vertex_t;
   vector<float>   *lostTrack_vertex_t;
   vector<float>   *pfCand_time;
   vector<float>   *lostTrack_time;
   vector<float>   *pfCand_timeError;
   vector<float>   *lostTrack_timeError;
   vector<int>     *pfCand_hasTrackDetails;
   vector<int>     *lostTrack_hasTrackDetails;
   vector<float>   *pfCand_dxy;
   vector<float>   *lostTrack_dxy;
   vector<float>   *pfCand_dxy_error;
   vector<float>   *lostTrack_dxy_error;
   vector<float>   *pfCand_dz;
   vector<float>   *lostTrack_dz;
   vector<float>   *pfCand_dz_error;
   vector<float>   *lostTrack_dz_error;
   vector<float>   *pfCand_track_pt;
   vector<float>   *lostTrack_track_pt;
   vector<float>   *pfCand_track_eta;
   vector<float>   *lostTrack_track_eta;
   vector<float>   *pfCand_track_phi;
   vector<float>   *lostTrack_track_phi;
   vector<float>   *pfCand_track_chi2;
   vector<float>   *lostTrack_track_chi2;
   vector<float>   *pfCand_track_ndof;
   vector<float>   *lostTrack_track_ndof;
   vector<float>   *pfCand_caloFraction;
   vector<float>   *lostTrack_caloFraction;
   vector<float>   *pfCand_hcalFraction;
   vector<float>   *lostTrack_hcalFraction;
   vector<float>   *pfCand_rawCaloFraction;
   vector<float>   *lostTrack_rawCaloFraction;
   vector<float>   *pfCand_rawHcalFraction;
   vector<float>   *lostTrack_rawHcalFraction;
   vector<int>     *ele_index;
   vector<float>   *ele_pt;
   vector<float>   *ele_eta;
   vector<float>   *ele_phi;
   vector<float>   *ele_mass;
   vector<float>   *ele_cc_ele_energy;
   vector<float>   *ele_cc_gamma_energy;
   vector<int>     *ele_cc_n_gamma;
   vector<float>   *ele_dxy;
   vector<float>   *ele_dxy_error;
   vector<float>   *ele_ip3d;
   vector<float>   *ele_trackMomentumAtVtx;
   vector<float>   *ele_trackMomentumAtCalo;
   vector<float>   *ele_trackMomentumOut;
   vector<float>   *ele_trackMomentumAtEleClus;
   vector<float>   *ele_trackMomentumAtVtxWithConstraint;
   vector<float>   *ele_ecalEnergy;
   vector<float>   *ele_ecalEnergy_error;
   vector<float>   *ele_eSuperClusterOverP;
   vector<float>   *ele_eSeedClusterOverP;
   vector<float>   *ele_eSeedClusterOverPout;
   vector<float>   *ele_eEleClusterOverPout;
   vector<float>   *ele_deltaEtaSuperClusterTrackAtVtx;
   vector<float>   *ele_deltaEtaSeedClusterTrackAtCalo;
   vector<float>   *ele_deltaEtaEleClusterTrackAtCalo;
   vector<float>   *ele_deltaEtaSeedClusterTrackAtVtx;
   vector<float>   *ele_deltaPhiEleClusterTrackAtCalo;
   vector<float>   *ele_deltaPhiSuperClusterTrackAtVtx;
   vector<float>   *ele_deltaPhiSeedClusterTrackAtCalo;
   vector<int>     *ele_mvaInput_earlyBrem;
   vector<int>     *ele_mvaInput_lateBrem;
   vector<float>   *ele_mvaInput_sigmaEtaEta;
   vector<float>   *ele_mvaInput_hadEnergy;
   vector<float>   *ele_mvaInput_deltaEta;
   vector<float>   *ele_gsfTrack_normalizedChi2;
   vector<int>     *ele_gsfTrack_numberOfValidHits;
   vector<float>   *ele_gsfTrack_pt;
   vector<float>   *ele_gsfTrack_pt_error;
   vector<float>   *ele_closestCtfTrack_normalizedChi2;
   vector<int>     *ele_closestCtfTrack_numberOfValidHits;
   vector<float>   *ele_sigmaEtaEta;
   vector<float>   *ele_sigmaIetaIeta;
   vector<float>   *ele_sigmaIphiIphi;
   vector<float>   *ele_sigmaIetaIphi;
   vector<float>   *ele_e1x5;
   vector<float>   *ele_e2x5Max;
   vector<float>   *ele_e5x5;
   vector<float>   *ele_r9;
   vector<float>   *ele_hcalDepth1OverEcal;
   vector<float>   *ele_hcalDepth2OverEcal;
   vector<float>   *ele_hcalDepth1OverEcalBc;
   vector<float>   *ele_hcalDepth2OverEcalBc;
   vector<float>   *ele_eLeft;
   vector<float>   *ele_eRight;
   vector<float>   *ele_eBottom;
   vector<float>   *ele_eTop;
   vector<float>   *ele_full5x5_sigmaEtaEta;
   vector<float>   *ele_full5x5_sigmaIetaIeta;
   vector<float>   *ele_full5x5_sigmaIphiIphi;
   vector<float>   *ele_full5x5_sigmaIetaIphi;
   vector<float>   *ele_full5x5_e1x5;
   vector<float>   *ele_full5x5_e2x5Max;
   vector<float>   *ele_full5x5_e5x5;
   vector<float>   *ele_full5x5_r9;
   vector<float>   *ele_full5x5_hcalDepth1OverEcal;
   vector<float>   *ele_full5x5_hcalDepth2OverEcal;
   vector<float>   *ele_full5x5_hcalDepth1OverEcalBc;
   vector<float>   *ele_full5x5_hcalDepth2OverEcalBc;
   vector<float>   *ele_full5x5_eLeft;
   vector<float>   *ele_full5x5_eRight;
   vector<float>   *ele_full5x5_eBottom;
   vector<float>   *ele_full5x5_eTop;
   vector<float>   *ele_full5x5_e2x5Left;
   vector<float>   *ele_full5x5_e2x5Right;
   vector<float>   *ele_full5x5_e2x5Bottom;
   vector<float>   *ele_full5x5_e2x5Top;
   vector<float>   *ele_hgcal_sigmaUU;
   vector<float>   *ele_hgcal_sigmaVV;
   vector<float>   *ele_hgcal_sigmaEE;
   vector<float>   *ele_hgcal_sigmaPP;
   vector<int>     *ele_hgcal_nLayers;
   vector<int>     *ele_hgcal_firstLayer;
   vector<int>     *ele_hgcal_lastLayer;
   vector<int>     *ele_hgcal_layerEfrac10;
   vector<int>     *ele_hgcal_layerEfrac90;
   vector<float>   *ele_hgcal_e4oEtot;
   vector<float>   *ele_hgcal_ecEnergy;
   vector<float>   *ele_hgcal_ecEnergyEE;
   vector<float>   *ele_hgcal_ecEnergyFH;
   vector<float>   *ele_hgcal_ecEnergyBH;
   vector<float>   *ele_hgcal_ecEt;
   vector<float>   *ele_hgcal_ecOrigEnergy;
   vector<float>   *ele_hgcal_ecOrigEt;
   vector<float>   *ele_hgcal_caloIsoRing0;
   vector<float>   *ele_hgcal_caloIsoRing1;
   vector<float>   *ele_hgcal_caloIsoRing2;
   vector<float>   *ele_hgcal_caloIsoRing3;
   vector<float>   *ele_hgcal_caloIsoRing4;
   vector<float>   *ele_hgcal_depthCompatibility;
   vector<float>   *ele_hgcal_expectedDepth;
   vector<float>   *ele_hgcal_expectedSigma;
   vector<float>   *ele_hgcal_measuredDepth;
   vector<float>   *ele_hgcal_pcaAxisX;
   vector<float>   *ele_hgcal_pcaAxisY;
   vector<float>   *ele_hgcal_pcaAxisZ;
   vector<float>   *ele_hgcal_pcaPositionX;
   vector<float>   *ele_hgcal_pcaPositionY;
   vector<float>   *ele_hgcal_pcaPositionZ;
   vector<float>   *ele_hgcal_pcaEig1;
   vector<float>   *ele_hgcal_pcaEig2;
   vector<float>   *ele_hgcal_pcaEig3;
   vector<float>   *ele_hgcal_pcaSig1;
   vector<float>   *ele_hgcal_pcaSig2;
   vector<float>   *ele_hgcal_pcaSig3;
   vector<int>     *muon_index;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_mass;
   vector<float>   *muon_dxy;
   vector<float>   *muon_dxy_error;
   vector<float>   *muon_normalizedChi2;
   vector<int>     *muon_numberOfValidHits;
   vector<float>   *muon_segmentCompatibility;
   vector<float>   *muon_caloCompatibility;
   vector<float>   *muon_pfEcalEnergy;
   vector<unsigned int> *muon_type;
   vector<int>     *muon_n_matches_DT_1;
   vector<int>     *muon_n_matches_DT_2;
   vector<int>     *muon_n_matches_DT_3;
   vector<int>     *muon_n_matches_DT_4;
   vector<int>     *muon_n_matches_CSC_1;
   vector<int>     *muon_n_matches_CSC_2;
   vector<int>     *muon_n_matches_CSC_3;
   vector<int>     *muon_n_matches_CSC_4;
   vector<int>     *muon_n_matches_RPC_1;
   vector<int>     *muon_n_matches_RPC_2;
   vector<int>     *muon_n_matches_RPC_3;
   vector<int>     *muon_n_matches_RPC_4;
   vector<int>     *muon_n_matches_GEM_1;
   vector<int>     *muon_n_matches_GEM_2;
   vector<int>     *muon_n_matches_GEM_3;
   vector<int>     *muon_n_matches_GEM_4;
   vector<int>     *muon_n_matches_ME0_1;
   vector<int>     *muon_n_matches_ME0_2;
   vector<int>     *muon_n_matches_ME0_3;
   vector<int>     *muon_n_matches_ME0_4;
   vector<int>     *muon_n_hits_DT_1;
   vector<int>     *muon_n_hits_DT_2;
   vector<int>     *muon_n_hits_DT_3;
   vector<int>     *muon_n_hits_DT_4;
   vector<int>     *muon_n_hits_CSC_1;
   vector<int>     *muon_n_hits_CSC_2;
   vector<int>     *muon_n_hits_CSC_3;
   vector<int>     *muon_n_hits_CSC_4;
   vector<int>     *muon_n_hits_RPC_1;
   vector<int>     *muon_n_hits_RPC_2;
   vector<int>     *muon_n_hits_RPC_3;
   vector<int>     *muon_n_hits_RPC_4;
   vector<int>     *muon_n_hits_GEM_1;
   vector<int>     *muon_n_hits_GEM_2;
   vector<int>     *muon_n_hits_GEM_3;
   vector<int>     *muon_n_hits_GEM_4;
   vector<int>     *muon_n_hits_ME0_1;
   vector<int>     *muon_n_hits_ME0_2;
   vector<int>     *muon_n_hits_ME0_3;
   vector<int>     *muon_n_hits_ME0_4;
   vector<int>     *isoTrack_index;
   vector<float>   *isoTrack_pt;
   vector<float>   *isoTrack_eta;
   vector<float>   *isoTrack_phi;
   vector<int>     *isoTrack_fromPV;
   vector<int>     *isoTrack_charge;
   vector<float>   *isoTrack_dxy;
   vector<float>   *isoTrack_dxy_error;
   vector<float>   *isoTrack_dz;
   vector<float>   *isoTrack_dz_error;
   vector<int>     *isoTrack_isHighPurityTrack;
   vector<int>     *isoTrack_isTightTrack;
   vector<int>     *isoTrack_isLooseTrack;
   vector<float>   *isoTrack_dEdxStrip;
   vector<float>   *isoTrack_dEdxPixel;
   vector<float>   *isoTrack_deltaEta;
   vector<float>   *isoTrack_deltaPhi;
   vector<int>     *isoTrack_n_ValidHits;
   vector<int>     *isoTrack_n_BadHits;
   vector<int>     *isoTrack_n_TimingHits;
   vector<int>     *isoTrack_n_ValidTimingHits;
   vector<int>     *isoTrack_n_LostTimingHits;
   vector<int>     *isoTrack_n_MuonHits;
   vector<int>     *isoTrack_n_ValidMuonHits;
   vector<int>     *isoTrack_n_LostMuonHits;
   vector<int>     *isoTrack_n_BadMuonHits;
   vector<int>     *isoTrack_n_ValidMuonDTHits;;
   vector<int>     *isoTrack_n_LostMuonDTHits;
   vector<int>     *isoTrack_n_BadMuonDTHits;
   vector<int>     *isoTrack_n_ValidMuonCSCHits;
   vector<int>     *isoTrack_n_LostMuonCSCHits;
   vector<int>     *isoTrack_n_BadMuonCSCHits;
   vector<int>     *isoTrack_n_ValidMuonRPCHits;
   vector<int>     *isoTrack_n_LostMuonRPCHits;
   vector<int>     *isoTrack_n_BadMuonRPCHits;
   vector<int>     *isoTrack_n_ValidMuonGEMHits;
   vector<int>     *isoTrack_n_LostMuonGEMHits;
   vector<int>     *isoTrack_n_BadMuonGEMHits;
   vector<int>     *isoTrack_n_ValidMuonME0Hits;
   vector<int>     *isoTrack_n_LostMuonME0Hits;
   vector<int>     *isoTrack_n_BadMuonME0Hits;
   vector<int>     *isoTrack_n_InactiveHits;
   vector<int>     *isoTrack_n_AllHits_TRACK;
   vector<int>     *isoTrack_n_AllHits_MISSING_INNER;
   vector<int>     *isoTrack_n_AllHits_MISSING_OUTER;
   vector<int>     *isoTrack_n_LostHits_TRACK;
   vector<int>     *isoTrack_n_LostHits_MISSING_INNER;
   vector<int>     *isoTrack_n_LostHits_MISSING_OUTER;
   vector<int>     *isoTrack_n_ValidPixelHits;
   vector<int>     *isoTrack_n_ValidStripHits;
   vector<int>     *isoTrack_n_LostPixelHits_TRACK;
   vector<int>     *isoTrack_n_LostPixelHits_MISSING_INNER;
   vector<int>     *isoTrack_n_LostPixelHits_MISSING_OUTER;
   vector<int>     *isoTrack_n_LostStripHits_TRACK;
   vector<int>     *isoTrack_n_LostStripHits_MISSING_INNER;
   vector<int>     *isoTrack_n_LostStripHits_MISSING_OUTER;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_genEventWeight;   //!
   TBranch        *b_trainingWeight;   //!
   TBranch        *b_sampleType;   //!
   TBranch        *b_tauType;   //!
   TBranch        *b_dataset_id;   //!
   TBranch        *b_dataset_group_id;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_t;   //!
   TBranch        *b_pv_xE;   //!
   TBranch        *b_pv_yE;   //!
   TBranch        *b_pv_zE;   //!
   TBranch        *b_pv_tE;   //!
   TBranch        *b_pv_chi2;   //!
   TBranch        *b_pv_ndof;   //!
   TBranch        *b_entry_index;   //!
   TBranch        *b_total_entries;   //!
   TBranch        *b_genLepton_index;   //!
   TBranch        *b_genLepton_kind;   //!
   TBranch        *b_genLepton_charge;   //!
   TBranch        *b_genLepton_vis_pt;   //!
   TBranch        *b_genLepton_vis_eta;   //!
   TBranch        *b_genLepton_vis_phi;   //!
   TBranch        *b_genLepton_vis_mass;   //!
   TBranch        *b_genLepton_lastMotherIndex;   //!
   TBranch        *b_genParticle_pdgId;   //!
   TBranch        *b_genParticle_mother;   //!
   TBranch        *b_genParticle_charge;   //!
   TBranch        *b_genParticle_isFirstCopy;   //!
   TBranch        *b_genParticle_isLastCopy;   //!
   TBranch        *b_genParticle_pt;   //!
   TBranch        *b_genParticle_eta;   //!
   TBranch        *b_genParticle_phi;   //!
   TBranch        *b_genParticle_mass;   //!
   TBranch        *b_genParticle_vtx_x;   //!
   TBranch        *b_genParticle_vtx_y;   //!
   TBranch        *b_genParticle_vtx_z;   //!
   TBranch        *b_genJet_index;   //!
   TBranch        *b_genJet_pt;   //!
   TBranch        *b_genJet_eta;   //!
   TBranch        *b_genJet_phi;   //!
   TBranch        *b_genJet_mass;   //!
   TBranch        *b_genJet_emEnergy;   //!
   TBranch        *b_genJet_hadEnergy;   //!
   TBranch        *b_genJet_invisibleEnergy;   //!
   TBranch        *b_genJet_auxiliaryEnergy;   //!
   TBranch        *b_genJet_chargedHadronEnergy;   //!
   TBranch        *b_genJet_neutralHadronEnergy;   //!
   TBranch        *b_genJet_chargedEmEnergy;   //!
   TBranch        *b_genJet_neutralEmEnergy;   //!
   TBranch        *b_genJet_muonEnergy;   //!
   TBranch        *b_genJet_chargedHadronMultiplicity;   //!
   TBranch        *b_genJet_neutralHadronMultiplicity;   //!
   TBranch        *b_genJet_chargedEmMultiplicity;   //!
   TBranch        *b_genJet_neutralEmMultiplicity;   //!
   TBranch        *b_genJet_muonMultiplicity;   //!
   TBranch        *b_genJet_n_bHadrons;   //!
   TBranch        *b_genJet_n_cHadrons;   //!
   TBranch        *b_genJet_n_partons;   //!
   TBranch        *b_genJet_n_leptons;   //!
   TBranch        *b_genJet_hadronFlavour;   //!
   TBranch        *b_genJet_partonFlavour;   //!
   TBranch        *b_jet_index;   //!
   TBranch        *b_fatJet_index;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_fatJet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_fatJet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_fatJet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_fatJet_mass;   //!
   TBranch        *b_jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_fatJet_neutralHadronEnergyFraction;   //!
   TBranch        *b_jet_neutralEmEnergyFraction;   //!
   TBranch        *b_fatJet_neutralEmEnergyFraction;   //!
   TBranch        *b_jet_nConstituents;   //!
   TBranch        *b_fatJet_nConstituents;   //!
   TBranch        *b_jet_chargedMultiplicity;   //!
   TBranch        *b_fatJet_chargedMultiplicity;   //!
   TBranch        *b_jet_neutralMultiplicity;   //!
   TBranch        *b_fatJet_neutralMultiplicity;   //!
   TBranch        *b_jet_partonFlavour;   //!
   TBranch        *b_fatJet_partonFlavour;   //!
   TBranch        *b_jet_hadronFlavour;   //!
   TBranch        *b_fatJet_hadronFlavour;   //!
   TBranch        *b_jet_m_softDrop;   //!
   TBranch        *b_fatJet_m_softDrop;   //!
   TBranch        *b_jet_nJettiness_tau1;   //!
   TBranch        *b_fatJet_nJettiness_tau1;   //!
   TBranch        *b_jet_nJettiness_tau2;   //!
   TBranch        *b_fatJet_nJettiness_tau2;   //!
   TBranch        *b_jet_nJettiness_tau3;   //!
   TBranch        *b_fatJet_nJettiness_tau3;   //!
   TBranch        *b_jet_nJettiness_tau4;   //!
   TBranch        *b_fatJet_nJettiness_tau4;   //!
   TBranch        *b_jet_subJet_pt;   //!
   TBranch        *b_fatJet_subJet_pt;   //!
   TBranch        *b_jet_subJet_eta;   //!
   TBranch        *b_fatJet_subJet_eta;   //!
   TBranch        *b_jet_subJet_phi;   //!
   TBranch        *b_fatJet_subJet_phi;   //!
   TBranch        *b_jet_subJet_mass;   //!
   TBranch        *b_fatJet_subJet_mass;   //!
   TBranch        *b_tau_index;   //!
   TBranch        *b_boostedTau_index;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_boostedTau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_boostedTau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_boostedTau_phi;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_boostedTau_mass;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_boostedTau_charge;   //!
   TBranch        *b_tau_decayMode;   //!
   TBranch        *b_boostedTau_decayMode;   //!
   TBranch        *b_tau_decayModeFinding;   //!
   TBranch        *b_boostedTau_decayModeFinding;   //!
   TBranch        *b_tau_decayModeFindingNewDMs;   //!
   TBranch        *b_boostedTau_decayModeFindingNewDMs;   //!
   TBranch        *b_tau_chargedIsoPtSum;   //!
   TBranch        *b_boostedTau_chargedIsoPtSum;   //!
   TBranch        *b_tau_chargedIsoPtSumdR03;   //!
   TBranch        *b_boostedTau_chargedIsoPtSumdR03;   //!
   TBranch        *b_tau_footprintCorrection;   //!
   TBranch        *b_boostedTau_footprintCorrection;   //!
   TBranch        *b_tau_footprintCorrectiondR03;   //!
   TBranch        *b_boostedTau_footprintCorrectiondR03;   //!
   TBranch        *b_tau_neutralIsoPtSum;   //!
   TBranch        *b_boostedTau_neutralIsoPtSum;   //!
   TBranch        *b_tau_neutralIsoPtSumWeight;   //!
   TBranch        *b_boostedTau_neutralIsoPtSumWeight;   //!
   TBranch        *b_tau_neutralIsoPtSumWeightdR03;   //!
   TBranch        *b_boostedTau_neutralIsoPtSumWeightdR03;   //!
   TBranch        *b_tau_neutralIsoPtSumdR03;   //!
   TBranch        *b_boostedTau_neutralIsoPtSumdR03;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_boostedTau_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_tau_photonPtSumOutsideSignalConedR03;   //!
   TBranch        *b_boostedTau_photonPtSumOutsideSignalConedR03;   //!
   TBranch        *b_tau_puCorrPtSum;   //!
   TBranch        *b_boostedTau_puCorrPtSum;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_boostedTau_byCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_tau_byCombinedIsolationDeltaBetaCorr3Hitsraw;   //!
   TBranch        *b_boostedTau_byCombinedIsolationDeltaBetaCorr3Hitsraw;   //!
   TBranch        *b_tau_byDeepTau2017v2p1VSe;   //!
   TBranch        *b_boostedTau_byDeepTau2017v2p1VSe;   //!
   TBranch        *b_tau_byDeepTau2017v2p1VSeraw;   //!
   TBranch        *b_boostedTau_byDeepTau2017v2p1VSeraw;   //!
   TBranch        *b_tau_byDeepTau2017v2p1VSmu;   //!
   TBranch        *b_boostedTau_byDeepTau2017v2p1VSmu;   //!
   TBranch        *b_tau_byDeepTau2017v2p1VSmuraw;   //!
   TBranch        *b_boostedTau_byDeepTau2017v2p1VSmuraw;   //!
   TBranch        *b_tau_byDeepTau2017v2p1VSjet;   //!
   TBranch        *b_boostedTau_byDeepTau2017v2p1VSjet;   //!
   TBranch        *b_tau_byDeepTau2017v2p1VSjetraw;   //!
   TBranch        *b_boostedTau_byDeepTau2017v2p1VSjetraw;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMwLT2017raw;   //!
   TBranch        *b_boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017raw;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBnewDMwLT2017raw;   //!
   TBranch        *b_boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017raw;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017;   //!
   TBranch        *b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw;   //!
   TBranch        *b_boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw;   //!
   TBranch        *b_tau_byIsolationMVADBnewDMwLTPhase2;   //!
   TBranch        *b_boostedTau_byIsolationMVADBnewDMwLTPhase2;   //!
   TBranch        *b_tau_byIsolationMVADBnewDMwLTPhase2raw;   //!
   TBranch        *b_boostedTau_byIsolationMVADBnewDMwLTPhase2raw;   //!
   TBranch        *b_tau_dxy_pca_x;   //!
   TBranch        *b_boostedTau_dxy_pca_x;   //!
   TBranch        *b_tau_dxy_pca_y;   //!
   TBranch        *b_boostedTau_dxy_pca_y;   //!
   TBranch        *b_tau_dxy_pca_z;   //!
   TBranch        *b_boostedTau_dxy_pca_z;   //!
   TBranch        *b_tau_dxy;   //!
   TBranch        *b_boostedTau_dxy;   //!
   TBranch        *b_tau_dxy_error;   //!
   TBranch        *b_boostedTau_dxy_error;   //!
   TBranch        *b_tau_ip3d;   //!
   TBranch        *b_boostedTau_ip3d;   //!
   TBranch        *b_tau_ip3d_error;   //!
   TBranch        *b_boostedTau_ip3d_error;   //!
   TBranch        *b_tau_dz;   //!
   TBranch        *b_boostedTau_dz;   //!
   TBranch        *b_tau_dz_error;   //!
   TBranch        *b_boostedTau_dz_error;   //!
   TBranch        *b_tau_hasSecondaryVertex;   //!
   TBranch        *b_boostedTau_hasSecondaryVertex;   //!
   TBranch        *b_tau_sv_x;   //!
   TBranch        *b_boostedTau_sv_x;   //!
   TBranch        *b_tau_sv_y;   //!
   TBranch        *b_boostedTau_sv_y;   //!
   TBranch        *b_tau_sv_z;   //!
   TBranch        *b_boostedTau_sv_z;   //!
   TBranch        *b_tau_flightLength_x;   //!
   TBranch        *b_boostedTau_flightLength_x;   //!
   TBranch        *b_tau_flightLength_y;   //!
   TBranch        *b_boostedTau_flightLength_y;   //!
   TBranch        *b_tau_flightLength_z;   //!
   TBranch        *b_boostedTau_flightLength_z;   //!
   TBranch        *b_tau_flightLength_sig;   //!
   TBranch        *b_boostedTau_flightLength_sig;   //!
   TBranch        *b_tau_pt_weighted_deta_strip;   //!
   TBranch        *b_boostedTau_pt_weighted_deta_strip;   //!
   TBranch        *b_tau_pt_weighted_dphi_strip;   //!
   TBranch        *b_boostedTau_pt_weighted_dphi_strip;   //!
   TBranch        *b_tau_pt_weighted_dr_signal;   //!
   TBranch        *b_boostedTau_pt_weighted_dr_signal;   //!
   TBranch        *b_tau_pt_weighted_dr_iso;   //!
   TBranch        *b_boostedTau_pt_weighted_dr_iso;   //!
   TBranch        *b_tau_leadingTrackNormChi2;   //!
   TBranch        *b_boostedTau_leadingTrackNormChi2;   //!
   TBranch        *b_tau_e_ratio;   //!
   TBranch        *b_boostedTau_e_ratio;   //!
   TBranch        *b_tau_gj_angle_diff;   //!
   TBranch        *b_boostedTau_gj_angle_diff;   //!
   TBranch        *b_tau_n_photons;   //!
   TBranch        *b_boostedTau_n_photons;   //!
   TBranch        *b_tau_emFraction;   //!
   TBranch        *b_boostedTau_emFraction;   //!
   TBranch        *b_tau_inside_ecal_crack;   //!
   TBranch        *b_boostedTau_inside_ecal_crack;   //!
   TBranch        *b_tau_leadChargedCand_etaAtEcalEntrance;   //!
   TBranch        *b_boostedTau_leadChargedCand_etaAtEcalEntrance;   //!
   TBranch        *b_pfCand_index;   //!
   TBranch        *b_lostTrack_index;   //!
   TBranch        *b_pfCand_tauSignal;   //!
   TBranch        *b_lostTrack_tauSignal;   //!
   TBranch        *b_pfCand_tauLeadChargedHadrCand;   //!
   TBranch        *b_lostTrack_tauLeadChargedHadrCand;   //!
   TBranch        *b_pfCand_tauIso;   //!
   TBranch        *b_lostTrack_tauIso;   //!
   TBranch        *b_pfCand_boostedTauSignal;   //!
   TBranch        *b_lostTrack_boostedTauSignal;   //!
   TBranch        *b_pfCand_boostedTauLeadChargedHadrCand;   //!
   TBranch        *b_lostTrack_boostedTauLeadChargedHadrCand;   //!
   TBranch        *b_pfCand_boostedTauIso;   //!
   TBranch        *b_lostTrack_boostedTauIso;   //!
   TBranch        *b_pfCand_jetDaughter;   //!
   TBranch        *b_lostTrack_jetDaughter;   //!
   TBranch        *b_pfCand_fatJetDaughter;   //!
   TBranch        *b_lostTrack_fatJetDaughter;   //!
   TBranch        *b_pfCand_subJetDaughter;   //!
   TBranch        *b_lostTrack_subJetDaughter;   //!
   TBranch        *b_pfCand_pt;   //!
   TBranch        *b_lostTrack_pt;   //!
   TBranch        *b_pfCand_eta;   //!
   TBranch        *b_lostTrack_eta;   //!
   TBranch        *b_pfCand_phi;   //!
   TBranch        *b_lostTrack_phi;   //!
   TBranch        *b_pfCand_mass;   //!
   TBranch        *b_lostTrack_mass;   //!
   TBranch        *b_pfCand_pvAssociationQuality;   //!
   TBranch        *b_lostTrack_pvAssociationQuality;   //!
   TBranch        *b_pfCand_fromPV;   //!
   TBranch        *b_lostTrack_fromPV;   //!
   TBranch        *b_pfCand_puppiWeight;   //!
   TBranch        *b_lostTrack_puppiWeight;   //!
   TBranch        *b_pfCand_puppiWeightNoLep;   //!
   TBranch        *b_lostTrack_puppiWeightNoLep;   //!
   TBranch        *b_pfCand_particleType;   //!
   TBranch        *b_lostTrack_particleType;   //!
   TBranch        *b_pfCand_charge;   //!
   TBranch        *b_lostTrack_charge;   //!
   TBranch        *b_pfCand_lostInnerHits;   //!
   TBranch        *b_lostTrack_lostInnerHits;   //!
   TBranch        *b_pfCand_nHits;   //!
   TBranch        *b_lostTrack_nHits;   //!
   TBranch        *b_pfCand_nPixelHits;   //!
   TBranch        *b_lostTrack_nPixelHits;   //!
   TBranch        *b_pfCand_nPixelLayers;   //!
   TBranch        *b_lostTrack_nPixelLayers;   //!
   TBranch        *b_pfCand_nStripLayers;   //!
   TBranch        *b_lostTrack_nStripLayers;   //!
   TBranch        *b_pfCand_vertex_x;   //!
   TBranch        *b_lostTrack_vertex_x;   //!
   TBranch        *b_pfCand_vertex_y;   //!
   TBranch        *b_lostTrack_vertex_y;   //!
   TBranch        *b_pfCand_vertex_z;   //!
   TBranch        *b_lostTrack_vertex_z;   //!
   TBranch        *b_pfCand_vertex_t;   //!
   TBranch        *b_lostTrack_vertex_t;   //!
   TBranch        *b_pfCand_time;   //!
   TBranch        *b_lostTrack_time;   //!
   TBranch        *b_pfCand_timeError;   //!
   TBranch        *b_lostTrack_timeError;   //!
   TBranch        *b_pfCand_hasTrackDetails;   //!
   TBranch        *b_lostTrack_hasTrackDetails;   //!
   TBranch        *b_pfCand_dxy;   //!
   TBranch        *b_lostTrack_dxy;   //!
   TBranch        *b_pfCand_dxy_error;   //!
   TBranch        *b_lostTrack_dxy_error;   //!
   TBranch        *b_pfCand_dz;   //!
   TBranch        *b_lostTrack_dz;   //!
   TBranch        *b_pfCand_dz_error;   //!
   TBranch        *b_lostTrack_dz_error;   //!
   TBranch        *b_pfCand_track_pt;   //!
   TBranch        *b_lostTrack_track_pt;   //!
   TBranch        *b_pfCand_track_eta;   //!
   TBranch        *b_lostTrack_track_eta;   //!
   TBranch        *b_pfCand_track_phi;   //!
   TBranch        *b_lostTrack_track_phi;   //!
   TBranch        *b_pfCand_track_chi2;   //!
   TBranch        *b_lostTrack_track_chi2;   //!
   TBranch        *b_pfCand_track_ndof;   //!
   TBranch        *b_lostTrack_track_ndof;   //!
   TBranch        *b_pfCand_caloFraction;   //!
   TBranch        *b_lostTrack_caloFraction;   //!
   TBranch        *b_pfCand_hcalFraction;   //!
   TBranch        *b_lostTrack_hcalFraction;   //!
   TBranch        *b_pfCand_rawCaloFraction;   //!
   TBranch        *b_lostTrack_rawCaloFraction;   //!
   TBranch        *b_pfCand_rawHcalFraction;   //!
   TBranch        *b_lostTrack_rawHcalFraction;   //!
   TBranch        *b_ele_index;   //!
   TBranch        *b_ele_pt;   //!
   TBranch        *b_ele_eta;   //!
   TBranch        *b_ele_phi;   //!
   TBranch        *b_ele_mass;   //!
   TBranch        *b_ele_cc_ele_energy;   //!
   TBranch        *b_ele_cc_gamma_energy;   //!
   TBranch        *b_ele_cc_n_gamma;   //!
   TBranch        *b_ele_dxy;   //!
   TBranch        *b_ele_dxy_error;   //!
   TBranch        *b_ele_ip3d;   //!
   TBranch        *b_ele_trackMomentumAtVtx;   //!
   TBranch        *b_ele_trackMomentumAtCalo;   //!
   TBranch        *b_ele_trackMomentumOut;   //!
   TBranch        *b_ele_trackMomentumAtEleClus;   //!
   TBranch        *b_ele_trackMomentumAtVtxWithConstraint;   //!
   TBranch        *b_ele_ecalEnergy;   //!
   TBranch        *b_ele_ecalEnergy_error;   //!
   TBranch        *b_ele_eSuperClusterOverP;   //!
   TBranch        *b_ele_eSeedClusterOverP;   //!
   TBranch        *b_ele_eSeedClusterOverPout;   //!
   TBranch        *b_ele_eEleClusterOverPout;   //!
   TBranch        *b_ele_deltaEtaSuperClusterTrackAtVtx;   //!
   TBranch        *b_ele_deltaEtaSeedClusterTrackAtCalo;   //!
   TBranch        *b_ele_deltaEtaEleClusterTrackAtCalo;   //!
   TBranch        *b_ele_deltaEtaSeedClusterTrackAtVtx;   //!
   TBranch        *b_ele_deltaPhiEleClusterTrackAtCalo;   //!
   TBranch        *b_ele_deltaPhiSuperClusterTrackAtVtx;   //!
   TBranch        *b_ele_deltaPhiSeedClusterTrackAtCalo;   //!
   TBranch        *b_ele_mvaInput_earlyBrem;   //!
   TBranch        *b_ele_mvaInput_lateBrem;   //!
   TBranch        *b_ele_mvaInput_sigmaEtaEta;   //!
   TBranch        *b_ele_mvaInput_hadEnergy;   //!
   TBranch        *b_ele_mvaInput_deltaEta;   //!
   TBranch        *b_ele_gsfTrack_normalizedChi2;   //!
   TBranch        *b_ele_gsfTrack_numberOfValidHits;   //!
   TBranch        *b_ele_gsfTrack_pt;   //!
   TBranch        *b_ele_gsfTrack_pt_error;   //!
   TBranch        *b_ele_closestCtfTrack_normalizedChi2;   //!
   TBranch        *b_ele_closestCtfTrack_numberOfValidHits;   //!
   TBranch        *b_ele_sigmaEtaEta;   //!
   TBranch        *b_ele_sigmaIetaIeta;   //!
   TBranch        *b_ele_sigmaIphiIphi;   //!
   TBranch        *b_ele_sigmaIetaIphi;   //!
   TBranch        *b_ele_e1x5;   //!
   TBranch        *b_ele_e2x5Max;   //!
   TBranch        *b_ele_e5x5;   //!
   TBranch        *b_ele_r9;   //!
   TBranch        *b_ele_hcalDepth1OverEcal;   //!
   TBranch        *b_ele_hcalDepth2OverEcal;   //!
   TBranch        *b_ele_hcalDepth1OverEcalBc;   //!
   TBranch        *b_ele_hcalDepth2OverEcalBc;   //!
   TBranch        *b_ele_eLeft;   //!
   TBranch        *b_ele_eRight;   //!
   TBranch        *b_ele_eBottom;   //!
   TBranch        *b_ele_eTop;   //!
   TBranch        *b_ele_full5x5_sigmaEtaEta;   //!
   TBranch        *b_ele_full5x5_sigmaIetaIeta;   //!
   TBranch        *b_ele_full5x5_sigmaIphiIphi;   //!
   TBranch        *b_ele_full5x5_sigmaIetaIphi;   //!
   TBranch        *b_ele_full5x5_e1x5;   //!
   TBranch        *b_ele_full5x5_e2x5Max;   //!
   TBranch        *b_ele_full5x5_e5x5;   //!
   TBranch        *b_ele_full5x5_r9;   //!
   TBranch        *b_ele_full5x5_hcalDepth1OverEcal;   //!
   TBranch        *b_ele_full5x5_hcalDepth2OverEcal;   //!
   TBranch        *b_ele_full5x5_hcalDepth1OverEcalBc;   //!
   TBranch        *b_ele_full5x5_hcalDepth2OverEcalBc;   //!
   TBranch        *b_ele_full5x5_eLeft;   //!
   TBranch        *b_ele_full5x5_eRight;   //!
   TBranch        *b_ele_full5x5_eBottom;   //!
   TBranch        *b_ele_full5x5_eTop;   //!
   TBranch        *b_ele_full5x5_e2x5Left;   //!
   TBranch        *b_ele_full5x5_e2x5Right;   //!
   TBranch        *b_ele_full5x5_e2x5Bottom;   //!
   TBranch        *b_ele_full5x5_e2x5Top;   //!
   TBranch        *b_ele_hgcal_sigmaUU;   //!
   TBranch        *b_ele_hgcal_sigmaVV;   //!
   TBranch        *b_ele_hgcal_sigmaEE;   //!
   TBranch        *b_ele_hgcal_sigmaPP;   //!
   TBranch        *b_ele_hgcal_nLayers;   //!
   TBranch        *b_ele_hgcal_firstLayer;   //!
   TBranch        *b_ele_hgcal_lastLayer;   //!
   TBranch        *b_ele_hgcal_layerEfrac10;   //!
   TBranch        *b_ele_hgcal_layerEfrac90;   //!
   TBranch        *b_ele_hgcal_e4oEtot;   //!
   TBranch        *b_ele_hgcal_ecEnergy;   //!
   TBranch        *b_ele_hgcal_ecEnergyEE;   //!
   TBranch        *b_ele_hgcal_ecEnergyFH;   //!
   TBranch        *b_ele_hgcal_ecEnergyBH;   //!
   TBranch        *b_ele_hgcal_ecEt;   //!
   TBranch        *b_ele_hgcal_ecOrigEnergy;   //!
   TBranch        *b_ele_hgcal_ecOrigEt;   //!
   TBranch        *b_ele_hgcal_caloIsoRing0;   //!
   TBranch        *b_ele_hgcal_caloIsoRing1;   //!
   TBranch        *b_ele_hgcal_caloIsoRing2;   //!
   TBranch        *b_ele_hgcal_caloIsoRing3;   //!
   TBranch        *b_ele_hgcal_caloIsoRing4;   //!
   TBranch        *b_ele_hgcal_depthCompatibility;   //!
   TBranch        *b_ele_hgcal_expectedDepth;   //!
   TBranch        *b_ele_hgcal_expectedSigma;   //!
   TBranch        *b_ele_hgcal_measuredDepth;   //!
   TBranch        *b_ele_hgcal_pcaAxisX;   //!
   TBranch        *b_ele_hgcal_pcaAxisY;   //!
   TBranch        *b_ele_hgcal_pcaAxisZ;   //!
   TBranch        *b_ele_hgcal_pcaPositionX;   //!
   TBranch        *b_ele_hgcal_pcaPositionY;   //!
   TBranch        *b_ele_hgcal_pcaPositionZ;   //!
   TBranch        *b_ele_hgcal_pcaEig1;   //!
   TBranch        *b_ele_hgcal_pcaEig2;   //!
   TBranch        *b_ele_hgcal_pcaEig3;   //!
   TBranch        *b_ele_hgcal_pcaSig1;   //!
   TBranch        *b_ele_hgcal_pcaSig2;   //!
   TBranch        *b_ele_hgcal_pcaSig3;   //!
   TBranch        *b_muon_index;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_mass;   //!
   TBranch        *b_muon_dxy;   //!
   TBranch        *b_muon_dxy_error;   //!
   TBranch        *b_muon_normalizedChi2;   //!
   TBranch        *b_muon_numberOfValidHits;   //!
   TBranch        *b_muon_segmentCompatibility;   //!
   TBranch        *b_muon_caloCompatibility;   //!
   TBranch        *b_muon_pfEcalEnergy;   //!
   TBranch        *b_muon_type;   //!
   TBranch        *b_muon_n_matches_DT_1;   //!
   TBranch        *b_muon_n_matches_DT_2;   //!
   TBranch        *b_muon_n_matches_DT_3;   //!
   TBranch        *b_muon_n_matches_DT_4;   //!
   TBranch        *b_muon_n_matches_CSC_1;   //!
   TBranch        *b_muon_n_matches_CSC_2;   //!
   TBranch        *b_muon_n_matches_CSC_3;   //!
   TBranch        *b_muon_n_matches_CSC_4;   //!
   TBranch        *b_muon_n_matches_RPC_1;   //!
   TBranch        *b_muon_n_matches_RPC_2;   //!
   TBranch        *b_muon_n_matches_RPC_3;   //!
   TBranch        *b_muon_n_matches_RPC_4;   //!
   TBranch        *b_muon_n_matches_GEM_1;   //!
   TBranch        *b_muon_n_matches_GEM_2;   //!
   TBranch        *b_muon_n_matches_GEM_3;   //!
   TBranch        *b_muon_n_matches_GEM_4;   //!
   TBranch        *b_muon_n_matches_ME0_1;   //!
   TBranch        *b_muon_n_matches_ME0_2;   //!
   TBranch        *b_muon_n_matches_ME0_3;   //!
   TBranch        *b_muon_n_matches_ME0_4;   //!
   TBranch        *b_muon_n_hits_DT_1;   //!
   TBranch        *b_muon_n_hits_DT_2;   //!
   TBranch        *b_muon_n_hits_DT_3;   //!
   TBranch        *b_muon_n_hits_DT_4;   //!
   TBranch        *b_muon_n_hits_CSC_1;   //!
   TBranch        *b_muon_n_hits_CSC_2;   //!
   TBranch        *b_muon_n_hits_CSC_3;   //!
   TBranch        *b_muon_n_hits_CSC_4;   //!
   TBranch        *b_muon_n_hits_RPC_1;   //!
   TBranch        *b_muon_n_hits_RPC_2;   //!
   TBranch        *b_muon_n_hits_RPC_3;   //!
   TBranch        *b_muon_n_hits_RPC_4;   //!
   TBranch        *b_muon_n_hits_GEM_1;   //!
   TBranch        *b_muon_n_hits_GEM_2;   //!
   TBranch        *b_muon_n_hits_GEM_3;   //!
   TBranch        *b_muon_n_hits_GEM_4;   //!
   TBranch        *b_muon_n_hits_ME0_1;   //!
   TBranch        *b_muon_n_hits_ME0_2;   //!
   TBranch        *b_muon_n_hits_ME0_3;   //!
   TBranch        *b_muon_n_hits_ME0_4;   //!
   TBranch        *b_isoTrack_index;   //!
   TBranch        *b_isoTrack_pt;   //!
   TBranch        *b_isoTrack_eta;   //!
   TBranch        *b_isoTrack_phi;   //!
   TBranch        *b_isoTrack_fromPV;   //!
   TBranch        *b_isoTrack_charge;   //!
   TBranch        *b_isoTrack_dxy;   //!
   TBranch        *b_isoTrack_dxy_error;   //!
   TBranch        *b_isoTrack_dz;   //!
   TBranch        *b_isoTrack_dz_error;   //!
   TBranch        *b_isoTrack_isHighPurityTrack;   //!
   TBranch        *b_isoTrack_isTightTrack;   //!
   TBranch        *b_isoTrack_isLooseTrack;   //!
   TBranch        *b_isoTrack_dEdxStrip;   //!
   TBranch        *b_isoTrack_dEdxPixel;   //!
   TBranch        *b_isoTrack_deltaEta;   //!
   TBranch        *b_isoTrack_deltaPhi;   //!
   TBranch        *b_isoTrack_n_ValidHits;   //!
   TBranch        *b_isoTrack_n_BadHits;   //!
   TBranch        *b_isoTrack_n_TimingHits;   //!
   TBranch        *b_isoTrack_n_ValidTimingHits;   //!
   TBranch        *b_isoTrack_n_LostTimingHits;   //!
   TBranch        *b_isoTrack_n_MuonHits;   //!
   TBranch        *b_isoTrack_n_ValidMuonHits;   //!
   TBranch        *b_isoTrack_n_LostMuonHits;   //!
   TBranch        *b_isoTrack_n_BadMuonHits;   //!
   TBranch        *b_isoTrack_n_ValidMuonDTHits;   //!
   TBranch        *b_isoTrack_n_LostMuonDTHits;   //!
   TBranch        *b_isoTrack_n_BadMuonDTHits;   //!
   TBranch        *b_isoTrack_n_ValidMuonCSCHits;   //!
   TBranch        *b_isoTrack_n_LostMuonCSCHits;   //!
   TBranch        *b_isoTrack_n_BadMuonCSCHits;   //!
   TBranch        *b_isoTrack_n_ValidMuonRPCHits;   //!
   TBranch        *b_isoTrack_n_LostMuonRPCHits;   //!
   TBranch        *b_isoTrack_n_BadMuonRPCHits;   //!
   TBranch        *b_isoTrack_n_ValidMuonGEMHits;   //!
   TBranch        *b_isoTrack_n_LostMuonGEMHits;   //!
   TBranch        *b_isoTrack_n_BadMuonGEMHits;   //!
   TBranch        *b_isoTrack_n_ValidMuonME0Hits;   //!
   TBranch        *b_isoTrack_n_LostMuonME0Hits;   //!
   TBranch        *b_isoTrack_n_BadMuonME0Hits;   //!
   TBranch        *b_isoTrack_n_InactiveHits;   //!
   TBranch        *b_isoTrack_n_AllHits_TRACK;   //!
   TBranch        *b_isoTrack_n_AllHits_MISSING_INNER;   //!
   TBranch        *b_isoTrack_n_AllHits_MISSING_OUTER;   //!
   TBranch        *b_isoTrack_n_LostHits_TRACK;   //!
   TBranch        *b_isoTrack_n_LostHits_MISSING_INNER;   //!
   TBranch        *b_isoTrack_n_LostHits_MISSING_OUTER;   //!
   TBranch        *b_isoTrack_n_ValidPixelHits;   //!
   TBranch        *b_isoTrack_n_ValidStripHits;   //!
   TBranch        *b_isoTrack_n_LostPixelHits_TRACK;   //!
   TBranch        *b_isoTrack_n_LostPixelHits_MISSING_INNER;   //!
   TBranch        *b_isoTrack_n_LostPixelHits_MISSING_OUTER;   //!
   TBranch        *b_isoTrack_n_LostStripHits_TRACK;   //!
   TBranch        *b_isoTrack_n_LostStripHits_MISSING_INNER;   //!
   TBranch        *b_isoTrack_n_LostStripHits_MISSING_OUTER;   //!

   MyTauClass(TTree *tree=0,TString filepath="/nfs/dust/cms/user/cardinia/TauReco/example_rootfiles/DYJetsToLL_M-50.root");
   virtual ~MyTauClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyTauClass_cxx
MyTauClass::MyTauClass(TTree *tree, TString filepath) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    if (filepath.Last('/')<0)
        inputName="Unknown";
    else
        inputName=filepath(filepath.Last('/')+1,filepath.Last('.')-filepath.Last('/')-1);
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filepath);
      if (!f || !f->IsOpen()) {
         f = new TFile(filepath);
      }
      f->GetObject("taus",tree);

   }
   Init(tree);
}

MyTauClass::~MyTauClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyTauClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyTauClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyTauClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genParticle_pdgId = 0;
   genParticle_mother = 0;
   genParticle_charge = 0;
   genParticle_isFirstCopy = 0;
   genParticle_isLastCopy = 0;
   genParticle_pt = 0;
   genParticle_eta = 0;
   genParticle_phi = 0;
   genParticle_mass = 0;
   genParticle_vtx_x = 0;
   genParticle_vtx_y = 0;
   genParticle_vtx_z = 0;
   jet_subJet_pt = 0;
   fatJet_subJet_pt = 0;
   jet_subJet_eta = 0;
   fatJet_subJet_eta = 0;
   jet_subJet_phi = 0;
   fatJet_subJet_phi = 0;
   jet_subJet_mass = 0;
   fatJet_subJet_mass = 0;
   pfCand_index = 0;
   lostTrack_index = 0;
   pfCand_tauSignal = 0;
   lostTrack_tauSignal = 0;
   pfCand_tauLeadChargedHadrCand = 0;
   lostTrack_tauLeadChargedHadrCand = 0;
   pfCand_tauIso = 0;
   lostTrack_tauIso = 0;
   pfCand_boostedTauSignal = 0;
   lostTrack_boostedTauSignal = 0;
   pfCand_boostedTauLeadChargedHadrCand = 0;
   lostTrack_boostedTauLeadChargedHadrCand = 0;
   pfCand_boostedTauIso = 0;
   lostTrack_boostedTauIso = 0;
   pfCand_jetDaughter = 0;
   lostTrack_jetDaughter = 0;
   pfCand_fatJetDaughter = 0;
   lostTrack_fatJetDaughter = 0;
   pfCand_subJetDaughter = 0;
   lostTrack_subJetDaughter = 0;
   pfCand_pt = 0;
   lostTrack_pt = 0;
   pfCand_eta = 0;
   lostTrack_eta = 0;
   pfCand_phi = 0;
   lostTrack_phi = 0;
   pfCand_mass = 0;
   lostTrack_mass = 0;
   pfCand_pvAssociationQuality = 0;
   lostTrack_pvAssociationQuality = 0;
   pfCand_fromPV = 0;
   lostTrack_fromPV = 0;
   pfCand_puppiWeight = 0;
   lostTrack_puppiWeight = 0;
   pfCand_puppiWeightNoLep = 0;
   lostTrack_puppiWeightNoLep = 0;
   pfCand_particleType = 0;
   lostTrack_particleType = 0;
   pfCand_charge = 0;
   lostTrack_charge = 0;
   pfCand_lostInnerHits = 0;
   lostTrack_lostInnerHits = 0;
   pfCand_nHits = 0;
   lostTrack_nHits = 0;
   pfCand_nPixelHits = 0;
   lostTrack_nPixelHits = 0;
   pfCand_nPixelLayers = 0;
   lostTrack_nPixelLayers = 0;
   pfCand_nStripLayers = 0;
   lostTrack_nStripLayers = 0;
   pfCand_vertex_x = 0;
   lostTrack_vertex_x = 0;
   pfCand_vertex_y = 0;
   lostTrack_vertex_y = 0;
   pfCand_vertex_z = 0;
   lostTrack_vertex_z = 0;
   pfCand_vertex_t = 0;
   lostTrack_vertex_t = 0;
   pfCand_time = 0;
   lostTrack_time = 0;
   pfCand_timeError = 0;
   lostTrack_timeError = 0;
   pfCand_hasTrackDetails = 0;
   lostTrack_hasTrackDetails = 0;
   pfCand_dxy = 0;
   lostTrack_dxy = 0;
   pfCand_dxy_error = 0;
   lostTrack_dxy_error = 0;
   pfCand_dz = 0;
   lostTrack_dz = 0;
   pfCand_dz_error = 0;
   lostTrack_dz_error = 0;
   pfCand_track_pt = 0;
   lostTrack_track_pt = 0;
   pfCand_track_eta = 0;
   lostTrack_track_eta = 0;
   pfCand_track_phi = 0;
   lostTrack_track_phi = 0;
   pfCand_track_chi2 = 0;
   lostTrack_track_chi2 = 0;
   pfCand_track_ndof = 0;
   lostTrack_track_ndof = 0;
   pfCand_caloFraction = 0;
   lostTrack_caloFraction = 0;
   pfCand_hcalFraction = 0;
   lostTrack_hcalFraction = 0;
   pfCand_rawCaloFraction = 0;
   lostTrack_rawCaloFraction = 0;
   pfCand_rawHcalFraction = 0;
   lostTrack_rawHcalFraction = 0;
   ele_index = 0;
   ele_pt = 0;
   ele_eta = 0;
   ele_phi = 0;
   ele_mass = 0;
   ele_cc_ele_energy = 0;
   ele_cc_gamma_energy = 0;
   ele_cc_n_gamma = 0;
   ele_dxy = 0;
   ele_dxy_error = 0;
   ele_ip3d = 0;
   ele_trackMomentumAtVtx = 0;
   ele_trackMomentumAtCalo = 0;
   ele_trackMomentumOut = 0;
   ele_trackMomentumAtEleClus = 0;
   ele_trackMomentumAtVtxWithConstraint = 0;
   ele_ecalEnergy = 0;
   ele_ecalEnergy_error = 0;
   ele_eSuperClusterOverP = 0;
   ele_eSeedClusterOverP = 0;
   ele_eSeedClusterOverPout = 0;
   ele_eEleClusterOverPout = 0;
   ele_deltaEtaSuperClusterTrackAtVtx = 0;
   ele_deltaEtaSeedClusterTrackAtCalo = 0;
   ele_deltaEtaEleClusterTrackAtCalo = 0;
   ele_deltaEtaSeedClusterTrackAtVtx = 0;
   ele_deltaPhiEleClusterTrackAtCalo = 0;
   ele_deltaPhiSuperClusterTrackAtVtx = 0;
   ele_deltaPhiSeedClusterTrackAtCalo = 0;
   ele_mvaInput_earlyBrem = 0;
   ele_mvaInput_lateBrem = 0;
   ele_mvaInput_sigmaEtaEta = 0;
   ele_mvaInput_hadEnergy = 0;
   ele_mvaInput_deltaEta = 0;
   ele_gsfTrack_normalizedChi2 = 0;
   ele_gsfTrack_numberOfValidHits = 0;
   ele_gsfTrack_pt = 0;
   ele_gsfTrack_pt_error = 0;
   ele_closestCtfTrack_normalizedChi2 = 0;
   ele_closestCtfTrack_numberOfValidHits = 0;
   ele_sigmaEtaEta = 0;
   ele_sigmaIetaIeta = 0;
   ele_sigmaIphiIphi = 0;
   ele_sigmaIetaIphi = 0;
   ele_e1x5 = 0;
   ele_e2x5Max = 0;
   ele_e5x5 = 0;
   ele_r9 = 0;
   ele_hcalDepth1OverEcal = 0;
   ele_hcalDepth2OverEcal = 0;
   ele_hcalDepth1OverEcalBc = 0;
   ele_hcalDepth2OverEcalBc = 0;
   ele_eLeft = 0;
   ele_eRight = 0;
   ele_eBottom = 0;
   ele_eTop = 0;
   ele_full5x5_sigmaEtaEta = 0;
   ele_full5x5_sigmaIetaIeta = 0;
   ele_full5x5_sigmaIphiIphi = 0;
   ele_full5x5_sigmaIetaIphi = 0;
   ele_full5x5_e1x5 = 0;
   ele_full5x5_e2x5Max = 0;
   ele_full5x5_e5x5 = 0;
   ele_full5x5_r9 = 0;
   ele_full5x5_hcalDepth1OverEcal = 0;
   ele_full5x5_hcalDepth2OverEcal = 0;
   ele_full5x5_hcalDepth1OverEcalBc = 0;
   ele_full5x5_hcalDepth2OverEcalBc = 0;
   ele_full5x5_eLeft = 0;
   ele_full5x5_eRight = 0;
   ele_full5x5_eBottom = 0;
   ele_full5x5_eTop = 0;
   ele_full5x5_e2x5Left = 0;
   ele_full5x5_e2x5Right = 0;
   ele_full5x5_e2x5Bottom = 0;
   ele_full5x5_e2x5Top = 0;
   ele_hgcal_sigmaUU = 0;
   ele_hgcal_sigmaVV = 0;
   ele_hgcal_sigmaEE = 0;
   ele_hgcal_sigmaPP = 0;
   ele_hgcal_nLayers = 0;
   ele_hgcal_firstLayer = 0;
   ele_hgcal_lastLayer = 0;
   ele_hgcal_layerEfrac10 = 0;
   ele_hgcal_layerEfrac90 = 0;
   ele_hgcal_e4oEtot = 0;
   ele_hgcal_ecEnergy = 0;
   ele_hgcal_ecEnergyEE = 0;
   ele_hgcal_ecEnergyFH = 0;
   ele_hgcal_ecEnergyBH = 0;
   ele_hgcal_ecEt = 0;
   ele_hgcal_ecOrigEnergy = 0;
   ele_hgcal_ecOrigEt = 0;
   ele_hgcal_caloIsoRing0 = 0;
   ele_hgcal_caloIsoRing1 = 0;
   ele_hgcal_caloIsoRing2 = 0;
   ele_hgcal_caloIsoRing3 = 0;
   ele_hgcal_caloIsoRing4 = 0;
   ele_hgcal_depthCompatibility = 0;
   ele_hgcal_expectedDepth = 0;
   ele_hgcal_expectedSigma = 0;
   ele_hgcal_measuredDepth = 0;
   ele_hgcal_pcaAxisX = 0;
   ele_hgcal_pcaAxisY = 0;
   ele_hgcal_pcaAxisZ = 0;
   ele_hgcal_pcaPositionX = 0;
   ele_hgcal_pcaPositionY = 0;
   ele_hgcal_pcaPositionZ = 0;
   ele_hgcal_pcaEig1 = 0;
   ele_hgcal_pcaEig2 = 0;
   ele_hgcal_pcaEig3 = 0;
   ele_hgcal_pcaSig1 = 0;
   ele_hgcal_pcaSig2 = 0;
   ele_hgcal_pcaSig3 = 0;
   muon_index = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_mass = 0;
   muon_dxy = 0;
   muon_dxy_error = 0;
   muon_normalizedChi2 = 0;
   muon_numberOfValidHits = 0;
   muon_segmentCompatibility = 0;
   muon_caloCompatibility = 0;
   muon_pfEcalEnergy = 0;
   muon_type = 0;
   muon_n_matches_DT_1 = 0;
   muon_n_matches_DT_2 = 0;
   muon_n_matches_DT_3 = 0;
   muon_n_matches_DT_4 = 0;
   muon_n_matches_CSC_1 = 0;
   muon_n_matches_CSC_2 = 0;
   muon_n_matches_CSC_3 = 0;
   muon_n_matches_CSC_4 = 0;
   muon_n_matches_RPC_1 = 0;
   muon_n_matches_RPC_2 = 0;
   muon_n_matches_RPC_3 = 0;
   muon_n_matches_RPC_4 = 0;
   muon_n_matches_GEM_1 = 0;
   muon_n_matches_GEM_2 = 0;
   muon_n_matches_GEM_3 = 0;
   muon_n_matches_GEM_4 = 0;
   muon_n_matches_ME0_1 = 0;
   muon_n_matches_ME0_2 = 0;
   muon_n_matches_ME0_3 = 0;
   muon_n_matches_ME0_4 = 0;
   muon_n_hits_DT_1 = 0;
   muon_n_hits_DT_2 = 0;
   muon_n_hits_DT_3 = 0;
   muon_n_hits_DT_4 = 0;
   muon_n_hits_CSC_1 = 0;
   muon_n_hits_CSC_2 = 0;
   muon_n_hits_CSC_3 = 0;
   muon_n_hits_CSC_4 = 0;
   muon_n_hits_RPC_1 = 0;
   muon_n_hits_RPC_2 = 0;
   muon_n_hits_RPC_3 = 0;
   muon_n_hits_RPC_4 = 0;
   muon_n_hits_GEM_1 = 0;
   muon_n_hits_GEM_2 = 0;
   muon_n_hits_GEM_3 = 0;
   muon_n_hits_GEM_4 = 0;
   muon_n_hits_ME0_1 = 0;
   muon_n_hits_ME0_2 = 0;
   muon_n_hits_ME0_3 = 0;
   muon_n_hits_ME0_4 = 0;
   isoTrack_index = 0;
   isoTrack_pt = 0;
   isoTrack_eta = 0;
   isoTrack_phi = 0;
   isoTrack_fromPV = 0;
   isoTrack_charge = 0;
   isoTrack_dxy = 0;
   isoTrack_dxy_error = 0;
   isoTrack_dz = 0;
   isoTrack_dz_error = 0;
   isoTrack_isHighPurityTrack = 0;
   isoTrack_isTightTrack = 0;
   isoTrack_isLooseTrack = 0;
   isoTrack_dEdxStrip = 0;
   isoTrack_dEdxPixel = 0;
   isoTrack_deltaEta = 0;
   isoTrack_deltaPhi = 0;
   isoTrack_n_ValidHits = 0;
   isoTrack_n_BadHits = 0;
   isoTrack_n_TimingHits = 0;
   isoTrack_n_ValidTimingHits = 0;
   isoTrack_n_LostTimingHits = 0;
   isoTrack_n_MuonHits = 0;
   isoTrack_n_ValidMuonHits = 0;
   isoTrack_n_LostMuonHits = 0;
   isoTrack_n_BadMuonHits = 0;
   isoTrack_n_ValidMuonDTHits = 0;
   isoTrack_n_LostMuonDTHits = 0;
   isoTrack_n_BadMuonDTHits = 0;
   isoTrack_n_ValidMuonCSCHits = 0;
   isoTrack_n_LostMuonCSCHits = 0;
   isoTrack_n_BadMuonCSCHits = 0;
   isoTrack_n_ValidMuonRPCHits = 0;
   isoTrack_n_LostMuonRPCHits = 0;
   isoTrack_n_BadMuonRPCHits = 0;
   isoTrack_n_ValidMuonGEMHits = 0;
   isoTrack_n_LostMuonGEMHits = 0;
   isoTrack_n_BadMuonGEMHits = 0;
   isoTrack_n_ValidMuonME0Hits = 0;
   isoTrack_n_LostMuonME0Hits = 0;
   isoTrack_n_BadMuonME0Hits = 0;
   isoTrack_n_InactiveHits = 0;
   isoTrack_n_AllHits_TRACK = 0;
   isoTrack_n_AllHits_MISSING_INNER = 0;
   isoTrack_n_AllHits_MISSING_OUTER = 0;
   isoTrack_n_LostHits_TRACK = 0;
   isoTrack_n_LostHits_MISSING_INNER = 0;
   isoTrack_n_LostHits_MISSING_OUTER = 0;
   isoTrack_n_ValidPixelHits = 0;
   isoTrack_n_ValidStripHits = 0;
   isoTrack_n_LostPixelHits_TRACK = 0;
   isoTrack_n_LostPixelHits_MISSING_INNER = 0;
   isoTrack_n_LostPixelHits_MISSING_OUTER = 0;
   isoTrack_n_LostStripHits_TRACK = 0;
   isoTrack_n_LostStripHits_MISSING_INNER = 0;
   isoTrack_n_LostStripHits_MISSING_OUTER = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
   fChain->SetBranchAddress("trainingWeight", &trainingWeight, &b_trainingWeight);
   fChain->SetBranchAddress("sampleType", &sampleType, &b_sampleType);
   fChain->SetBranchAddress("tauType", &tauType, &b_tauType);
   fChain->SetBranchAddress("dataset_id", &dataset_id, &b_dataset_id);
   fChain->SetBranchAddress("dataset_group_id", &dataset_group_id, &b_dataset_group_id);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
   fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_t", &pv_t, &b_pv_t);
   fChain->SetBranchAddress("pv_xE", &pv_xE, &b_pv_xE);
   fChain->SetBranchAddress("pv_yE", &pv_yE, &b_pv_yE);
   fChain->SetBranchAddress("pv_zE", &pv_zE, &b_pv_zE);
   fChain->SetBranchAddress("pv_tE", &pv_tE, &b_pv_tE);
   fChain->SetBranchAddress("pv_chi2", &pv_chi2, &b_pv_chi2);
   fChain->SetBranchAddress("pv_ndof", &pv_ndof, &b_pv_ndof);
   fChain->SetBranchAddress("entry_index", &entry_index, &b_entry_index);
   fChain->SetBranchAddress("total_entries", &total_entries, &b_total_entries);
   fChain->SetBranchAddress("genLepton_index", &genLepton_index, &b_genLepton_index);
   fChain->SetBranchAddress("genLepton_kind", &genLepton_kind, &b_genLepton_kind);
   fChain->SetBranchAddress("genLepton_charge", &genLepton_charge, &b_genLepton_charge);
   fChain->SetBranchAddress("genLepton_vis_pt", &genLepton_vis_pt, &b_genLepton_vis_pt);
   fChain->SetBranchAddress("genLepton_vis_eta", &genLepton_vis_eta, &b_genLepton_vis_eta);
   fChain->SetBranchAddress("genLepton_vis_phi", &genLepton_vis_phi, &b_genLepton_vis_phi);
   fChain->SetBranchAddress("genLepton_vis_mass", &genLepton_vis_mass, &b_genLepton_vis_mass);
   fChain->SetBranchAddress("genLepton_lastMotherIndex", &genLepton_lastMotherIndex, &b_genLepton_lastMotherIndex);
   fChain->SetBranchAddress("genParticle_pdgId", &genParticle_pdgId, &b_genParticle_pdgId);
   fChain->SetBranchAddress("genParticle_mother", &genParticle_mother, &b_genParticle_mother);
   fChain->SetBranchAddress("genParticle_charge", &genParticle_charge, &b_genParticle_charge);
   fChain->SetBranchAddress("genParticle_isFirstCopy", &genParticle_isFirstCopy, &b_genParticle_isFirstCopy);
   fChain->SetBranchAddress("genParticle_isLastCopy", &genParticle_isLastCopy, &b_genParticle_isLastCopy);
   fChain->SetBranchAddress("genParticle_pt", &genParticle_pt, &b_genParticle_pt);
   fChain->SetBranchAddress("genParticle_eta", &genParticle_eta, &b_genParticle_eta);
   fChain->SetBranchAddress("genParticle_phi", &genParticle_phi, &b_genParticle_phi);
   fChain->SetBranchAddress("genParticle_mass", &genParticle_mass, &b_genParticle_mass);
   fChain->SetBranchAddress("genParticle_vtx_x", &genParticle_vtx_x, &b_genParticle_vtx_x);
   fChain->SetBranchAddress("genParticle_vtx_y", &genParticle_vtx_y, &b_genParticle_vtx_y);
   fChain->SetBranchAddress("genParticle_vtx_z", &genParticle_vtx_z, &b_genParticle_vtx_z);
   fChain->SetBranchAddress("genJet_index", &genJet_index, &b_genJet_index);
   fChain->SetBranchAddress("genJet_pt", &genJet_pt, &b_genJet_pt);
   fChain->SetBranchAddress("genJet_eta", &genJet_eta, &b_genJet_eta);
   fChain->SetBranchAddress("genJet_phi", &genJet_phi, &b_genJet_phi);
   fChain->SetBranchAddress("genJet_mass", &genJet_mass, &b_genJet_mass);
   fChain->SetBranchAddress("genJet_emEnergy", &genJet_emEnergy, &b_genJet_emEnergy);
   fChain->SetBranchAddress("genJet_hadEnergy", &genJet_hadEnergy, &b_genJet_hadEnergy);
   fChain->SetBranchAddress("genJet_invisibleEnergy", &genJet_invisibleEnergy, &b_genJet_invisibleEnergy);
   fChain->SetBranchAddress("genJet_auxiliaryEnergy", &genJet_auxiliaryEnergy, &b_genJet_auxiliaryEnergy);
   fChain->SetBranchAddress("genJet_chargedHadronEnergy", &genJet_chargedHadronEnergy, &b_genJet_chargedHadronEnergy);
   fChain->SetBranchAddress("genJet_neutralHadronEnergy", &genJet_neutralHadronEnergy, &b_genJet_neutralHadronEnergy);
   fChain->SetBranchAddress("genJet_chargedEmEnergy", &genJet_chargedEmEnergy, &b_genJet_chargedEmEnergy);
   fChain->SetBranchAddress("genJet_neutralEmEnergy", &genJet_neutralEmEnergy, &b_genJet_neutralEmEnergy);
   fChain->SetBranchAddress("genJet_muonEnergy", &genJet_muonEnergy, &b_genJet_muonEnergy);
   fChain->SetBranchAddress("genJet_chargedHadronMultiplicity", &genJet_chargedHadronMultiplicity, &b_genJet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("genJet_neutralHadronMultiplicity", &genJet_neutralHadronMultiplicity, &b_genJet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("genJet_chargedEmMultiplicity", &genJet_chargedEmMultiplicity, &b_genJet_chargedEmMultiplicity);
   fChain->SetBranchAddress("genJet_neutralEmMultiplicity", &genJet_neutralEmMultiplicity, &b_genJet_neutralEmMultiplicity);
   fChain->SetBranchAddress("genJet_muonMultiplicity", &genJet_muonMultiplicity, &b_genJet_muonMultiplicity);
   fChain->SetBranchAddress("genJet_n_bHadrons", &genJet_n_bHadrons, &b_genJet_n_bHadrons);
   fChain->SetBranchAddress("genJet_n_cHadrons", &genJet_n_cHadrons, &b_genJet_n_cHadrons);
   fChain->SetBranchAddress("genJet_n_partons", &genJet_n_partons, &b_genJet_n_partons);
   fChain->SetBranchAddress("genJet_n_leptons", &genJet_n_leptons, &b_genJet_n_leptons);
   fChain->SetBranchAddress("genJet_hadronFlavour", &genJet_hadronFlavour, &b_genJet_hadronFlavour);
   fChain->SetBranchAddress("genJet_partonFlavour", &genJet_partonFlavour, &b_genJet_partonFlavour);
   fChain->SetBranchAddress("jet_index", &jet_index, &b_jet_index);
   fChain->SetBranchAddress("fatJet_index", &fatJet_index, &b_fatJet_index);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("fatJet_pt", &fatJet_pt, &b_fatJet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("fatJet_eta", &fatJet_eta, &b_fatJet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("fatJet_phi", &fatJet_phi, &b_fatJet_phi);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("fatJet_mass", &fatJet_mass, &b_fatJet_mass);
   fChain->SetBranchAddress("jet_neutralHadronEnergyFraction", &jet_neutralHadronEnergyFraction, &b_jet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("fatJet_neutralHadronEnergyFraction", &fatJet_neutralHadronEnergyFraction, &b_fatJet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("jet_neutralEmEnergyFraction", &jet_neutralEmEnergyFraction, &b_jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("fatJet_neutralEmEnergyFraction", &fatJet_neutralEmEnergyFraction, &b_fatJet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("jet_nConstituents", &jet_nConstituents, &b_jet_nConstituents);
   fChain->SetBranchAddress("fatJet_nConstituents", &fatJet_nConstituents, &b_fatJet_nConstituents);
   fChain->SetBranchAddress("jet_chargedMultiplicity", &jet_chargedMultiplicity, &b_jet_chargedMultiplicity);
   fChain->SetBranchAddress("fatJet_chargedMultiplicity", &fatJet_chargedMultiplicity, &b_fatJet_chargedMultiplicity);
   fChain->SetBranchAddress("jet_neutralMultiplicity", &jet_neutralMultiplicity, &b_jet_neutralMultiplicity);
   fChain->SetBranchAddress("fatJet_neutralMultiplicity", &fatJet_neutralMultiplicity, &b_fatJet_neutralMultiplicity);
   fChain->SetBranchAddress("jet_partonFlavour", &jet_partonFlavour, &b_jet_partonFlavour);
   fChain->SetBranchAddress("fatJet_partonFlavour", &fatJet_partonFlavour, &b_fatJet_partonFlavour);
   fChain->SetBranchAddress("jet_hadronFlavour", &jet_hadronFlavour, &b_jet_hadronFlavour);
   fChain->SetBranchAddress("fatJet_hadronFlavour", &fatJet_hadronFlavour, &b_fatJet_hadronFlavour);
   fChain->SetBranchAddress("jet_m_softDrop", &jet_m_softDrop, &b_jet_m_softDrop);
   fChain->SetBranchAddress("fatJet_m_softDrop", &fatJet_m_softDrop, &b_fatJet_m_softDrop);
   fChain->SetBranchAddress("jet_nJettiness_tau1", &jet_nJettiness_tau1, &b_jet_nJettiness_tau1);
   fChain->SetBranchAddress("fatJet_nJettiness_tau1", &fatJet_nJettiness_tau1, &b_fatJet_nJettiness_tau1);
   fChain->SetBranchAddress("jet_nJettiness_tau2", &jet_nJettiness_tau2, &b_jet_nJettiness_tau2);
   fChain->SetBranchAddress("fatJet_nJettiness_tau2", &fatJet_nJettiness_tau2, &b_fatJet_nJettiness_tau2);
   fChain->SetBranchAddress("jet_nJettiness_tau3", &jet_nJettiness_tau3, &b_jet_nJettiness_tau3);
   fChain->SetBranchAddress("fatJet_nJettiness_tau3", &fatJet_nJettiness_tau3, &b_fatJet_nJettiness_tau3);
   fChain->SetBranchAddress("jet_nJettiness_tau4", &jet_nJettiness_tau4, &b_jet_nJettiness_tau4);
   fChain->SetBranchAddress("fatJet_nJettiness_tau4", &fatJet_nJettiness_tau4, &b_fatJet_nJettiness_tau4);
   fChain->SetBranchAddress("jet_subJet_pt", &jet_subJet_pt, &b_jet_subJet_pt);
   fChain->SetBranchAddress("fatJet_subJet_pt", &fatJet_subJet_pt, &b_fatJet_subJet_pt);
   fChain->SetBranchAddress("jet_subJet_eta", &jet_subJet_eta, &b_jet_subJet_eta);
   fChain->SetBranchAddress("fatJet_subJet_eta", &fatJet_subJet_eta, &b_fatJet_subJet_eta);
   fChain->SetBranchAddress("jet_subJet_phi", &jet_subJet_phi, &b_jet_subJet_phi);
   fChain->SetBranchAddress("fatJet_subJet_phi", &fatJet_subJet_phi, &b_fatJet_subJet_phi);
   fChain->SetBranchAddress("jet_subJet_mass", &jet_subJet_mass, &b_jet_subJet_mass);
   fChain->SetBranchAddress("fatJet_subJet_mass", &fatJet_subJet_mass, &b_fatJet_subJet_mass);
   fChain->SetBranchAddress("tau_index", &tau_index, &b_tau_index);
   fChain->SetBranchAddress("boostedTau_index", &boostedTau_index, &b_boostedTau_index);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("boostedTau_pt", &boostedTau_pt, &b_boostedTau_pt);
   fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("boostedTau_eta", &boostedTau_eta, &b_boostedTau_eta);
   fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("boostedTau_phi", &boostedTau_phi, &b_boostedTau_phi);
   fChain->SetBranchAddress("tau_mass", &tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("boostedTau_mass", &boostedTau_mass, &b_boostedTau_mass);
   fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("boostedTau_charge", &boostedTau_charge, &b_boostedTau_charge);
   fChain->SetBranchAddress("tau_decayMode", &tau_decayMode, &b_tau_decayMode);
   fChain->SetBranchAddress("boostedTau_decayMode", &boostedTau_decayMode, &b_boostedTau_decayMode);
   fChain->SetBranchAddress("tau_decayModeFinding", &tau_decayModeFinding, &b_tau_decayModeFinding);
   fChain->SetBranchAddress("boostedTau_decayModeFinding", &boostedTau_decayModeFinding, &b_boostedTau_decayModeFinding);
   fChain->SetBranchAddress("tau_decayModeFindingNewDMs", &tau_decayModeFindingNewDMs, &b_tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("boostedTau_decayModeFindingNewDMs", &boostedTau_decayModeFindingNewDMs, &b_boostedTau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("tau_chargedIsoPtSum", &tau_chargedIsoPtSum, &b_tau_chargedIsoPtSum);
   fChain->SetBranchAddress("boostedTau_chargedIsoPtSum", &boostedTau_chargedIsoPtSum, &b_boostedTau_chargedIsoPtSum);
   fChain->SetBranchAddress("tau_chargedIsoPtSumdR03", &tau_chargedIsoPtSumdR03, &b_tau_chargedIsoPtSumdR03);
   fChain->SetBranchAddress("boostedTau_chargedIsoPtSumdR03", &boostedTau_chargedIsoPtSumdR03, &b_boostedTau_chargedIsoPtSumdR03);
   fChain->SetBranchAddress("tau_footprintCorrection", &tau_footprintCorrection, &b_tau_footprintCorrection);
   fChain->SetBranchAddress("boostedTau_footprintCorrection", &boostedTau_footprintCorrection, &b_boostedTau_footprintCorrection);
   fChain->SetBranchAddress("tau_footprintCorrectiondR03", &tau_footprintCorrectiondR03, &b_tau_footprintCorrectiondR03);
   fChain->SetBranchAddress("boostedTau_footprintCorrectiondR03", &boostedTau_footprintCorrectiondR03, &b_boostedTau_footprintCorrectiondR03);
   fChain->SetBranchAddress("tau_neutralIsoPtSum", &tau_neutralIsoPtSum, &b_tau_neutralIsoPtSum);
   fChain->SetBranchAddress("boostedTau_neutralIsoPtSum", &boostedTau_neutralIsoPtSum, &b_boostedTau_neutralIsoPtSum);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeight", &tau_neutralIsoPtSumWeight, &b_tau_neutralIsoPtSumWeight);
   fChain->SetBranchAddress("boostedTau_neutralIsoPtSumWeight", &boostedTau_neutralIsoPtSumWeight, &b_boostedTau_neutralIsoPtSumWeight);
   fChain->SetBranchAddress("tau_neutralIsoPtSumWeightdR03", &tau_neutralIsoPtSumWeightdR03, &b_tau_neutralIsoPtSumWeightdR03);
   fChain->SetBranchAddress("boostedTau_neutralIsoPtSumWeightdR03", &boostedTau_neutralIsoPtSumWeightdR03, &b_boostedTau_neutralIsoPtSumWeightdR03);
   fChain->SetBranchAddress("tau_neutralIsoPtSumdR03", &tau_neutralIsoPtSumdR03, &b_tau_neutralIsoPtSumdR03);
   fChain->SetBranchAddress("boostedTau_neutralIsoPtSumdR03", &boostedTau_neutralIsoPtSumdR03, &b_boostedTau_neutralIsoPtSumdR03);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalCone", &tau_photonPtSumOutsideSignalCone, &b_tau_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("boostedTau_photonPtSumOutsideSignalCone", &boostedTau_photonPtSumOutsideSignalCone, &b_boostedTau_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("tau_photonPtSumOutsideSignalConedR03", &tau_photonPtSumOutsideSignalConedR03, &b_tau_photonPtSumOutsideSignalConedR03);
   fChain->SetBranchAddress("boostedTau_photonPtSumOutsideSignalConedR03", &boostedTau_photonPtSumOutsideSignalConedR03, &b_boostedTau_photonPtSumOutsideSignalConedR03);
   fChain->SetBranchAddress("tau_puCorrPtSum", &tau_puCorrPtSum, &b_tau_puCorrPtSum);
   fChain->SetBranchAddress("boostedTau_puCorrPtSum", &boostedTau_puCorrPtSum, &b_boostedTau_puCorrPtSum);
   fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorr3Hits", &tau_byCombinedIsolationDeltaBetaCorr3Hits, &b_tau_byCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("boostedTau_byCombinedIsolationDeltaBetaCorr3Hits", &boostedTau_byCombinedIsolationDeltaBetaCorr3Hits, &b_boostedTau_byCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("tau_byCombinedIsolationDeltaBetaCorr3Hitsraw", &tau_byCombinedIsolationDeltaBetaCorr3Hitsraw, &b_tau_byCombinedIsolationDeltaBetaCorr3Hitsraw);
   fChain->SetBranchAddress("boostedTau_byCombinedIsolationDeltaBetaCorr3Hitsraw", &boostedTau_byCombinedIsolationDeltaBetaCorr3Hitsraw, &b_boostedTau_byCombinedIsolationDeltaBetaCorr3Hitsraw);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSe", &tau_byDeepTau2017v2p1VSe, &b_tau_byDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("boostedTau_byDeepTau2017v2p1VSe", &boostedTau_byDeepTau2017v2p1VSe, &b_boostedTau_byDeepTau2017v2p1VSe);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSeraw", &tau_byDeepTau2017v2p1VSeraw, &b_tau_byDeepTau2017v2p1VSeraw);
   fChain->SetBranchAddress("boostedTau_byDeepTau2017v2p1VSeraw", &boostedTau_byDeepTau2017v2p1VSeraw, &b_boostedTau_byDeepTau2017v2p1VSeraw);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSmu", &tau_byDeepTau2017v2p1VSmu, &b_tau_byDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("boostedTau_byDeepTau2017v2p1VSmu", &boostedTau_byDeepTau2017v2p1VSmu, &b_boostedTau_byDeepTau2017v2p1VSmu);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSmuraw", &tau_byDeepTau2017v2p1VSmuraw, &b_tau_byDeepTau2017v2p1VSmuraw);
   fChain->SetBranchAddress("boostedTau_byDeepTau2017v2p1VSmuraw", &boostedTau_byDeepTau2017v2p1VSmuraw, &b_boostedTau_byDeepTau2017v2p1VSmuraw);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSjet", &tau_byDeepTau2017v2p1VSjet, &b_tau_byDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("boostedTau_byDeepTau2017v2p1VSjet", &boostedTau_byDeepTau2017v2p1VSjet, &b_boostedTau_byDeepTau2017v2p1VSjet);
   fChain->SetBranchAddress("tau_byDeepTau2017v2p1VSjetraw", &tau_byDeepTau2017v2p1VSjetraw, &b_tau_byDeepTau2017v2p1VSjetraw);
   fChain->SetBranchAddress("boostedTau_byDeepTau2017v2p1VSjetraw", &boostedTau_byDeepTau2017v2p1VSjetraw, &b_boostedTau_byDeepTau2017v2p1VSjetraw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLT2017", &tau_byIsolationMVArun2017v2DBoldDMwLT2017, &b_tau_byIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017", &boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017, &b_boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMwLT2017raw", &tau_byIsolationMVArun2017v2DBoldDMwLT2017raw, &b_tau_byIsolationMVArun2017v2DBoldDMwLT2017raw);
   fChain->SetBranchAddress("boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017raw", &boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017raw, &b_boostedTau_byIsolationMVArun2017v2DBoldDMwLT2017raw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBnewDMwLT2017", &tau_byIsolationMVArun2017v2DBnewDMwLT2017, &b_tau_byIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017", &boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017, &b_boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBnewDMwLT2017raw", &tau_byIsolationMVArun2017v2DBnewDMwLT2017raw, &b_tau_byIsolationMVArun2017v2DBnewDMwLT2017raw);
   fChain->SetBranchAddress("boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017raw", &boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017raw, &b_boostedTau_byIsolationMVArun2017v2DBnewDMwLT2017raw);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017", &boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017, &b_boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017);
   fChain->SetBranchAddress("tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw", &tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw, &b_tau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw);
   fChain->SetBranchAddress("boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw", &boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw, &b_boostedTau_byIsolationMVArun2017v2DBoldDMdR0p3wLT2017raw);
   fChain->SetBranchAddress("tau_byIsolationMVADBnewDMwLTPhase2", &tau_byIsolationMVADBnewDMwLTPhase2, &b_tau_byIsolationMVADBnewDMwLTPhase2);
   fChain->SetBranchAddress("boostedTau_byIsolationMVADBnewDMwLTPhase2", &boostedTau_byIsolationMVADBnewDMwLTPhase2, &b_boostedTau_byIsolationMVADBnewDMwLTPhase2);
   fChain->SetBranchAddress("tau_byIsolationMVADBnewDMwLTPhase2raw", &tau_byIsolationMVADBnewDMwLTPhase2raw, &b_tau_byIsolationMVADBnewDMwLTPhase2raw);
   fChain->SetBranchAddress("boostedTau_byIsolationMVADBnewDMwLTPhase2raw", &boostedTau_byIsolationMVADBnewDMwLTPhase2raw, &b_boostedTau_byIsolationMVADBnewDMwLTPhase2raw);
   fChain->SetBranchAddress("tau_dxy_pca_x", &tau_dxy_pca_x, &b_tau_dxy_pca_x);
   fChain->SetBranchAddress("boostedTau_dxy_pca_x", &boostedTau_dxy_pca_x, &b_boostedTau_dxy_pca_x);
   fChain->SetBranchAddress("tau_dxy_pca_y", &tau_dxy_pca_y, &b_tau_dxy_pca_y);
   fChain->SetBranchAddress("boostedTau_dxy_pca_y", &boostedTau_dxy_pca_y, &b_boostedTau_dxy_pca_y);
   fChain->SetBranchAddress("tau_dxy_pca_z", &tau_dxy_pca_z, &b_tau_dxy_pca_z);
   fChain->SetBranchAddress("boostedTau_dxy_pca_z", &boostedTau_dxy_pca_z, &b_boostedTau_dxy_pca_z);
   fChain->SetBranchAddress("tau_dxy", &tau_dxy, &b_tau_dxy);
   fChain->SetBranchAddress("boostedTau_dxy", &boostedTau_dxy, &b_boostedTau_dxy);
   fChain->SetBranchAddress("tau_dxy_error", &tau_dxy_error, &b_tau_dxy_error);
   fChain->SetBranchAddress("boostedTau_dxy_error", &boostedTau_dxy_error, &b_boostedTau_dxy_error);
   fChain->SetBranchAddress("tau_ip3d", &tau_ip3d, &b_tau_ip3d);
   fChain->SetBranchAddress("boostedTau_ip3d", &boostedTau_ip3d, &b_boostedTau_ip3d);
   fChain->SetBranchAddress("tau_ip3d_error", &tau_ip3d_error, &b_tau_ip3d_error);
   fChain->SetBranchAddress("boostedTau_ip3d_error", &boostedTau_ip3d_error, &b_boostedTau_ip3d_error);
   fChain->SetBranchAddress("tau_dz", &tau_dz, &b_tau_dz);
   fChain->SetBranchAddress("boostedTau_dz", &boostedTau_dz, &b_boostedTau_dz);
   fChain->SetBranchAddress("tau_dz_error", &tau_dz_error, &b_tau_dz_error);
   fChain->SetBranchAddress("boostedTau_dz_error", &boostedTau_dz_error, &b_boostedTau_dz_error);
   fChain->SetBranchAddress("tau_hasSecondaryVertex", &tau_hasSecondaryVertex, &b_tau_hasSecondaryVertex);
   fChain->SetBranchAddress("boostedTau_hasSecondaryVertex", &boostedTau_hasSecondaryVertex, &b_boostedTau_hasSecondaryVertex);
   fChain->SetBranchAddress("tau_sv_x", &tau_sv_x, &b_tau_sv_x);
   fChain->SetBranchAddress("boostedTau_sv_x", &boostedTau_sv_x, &b_boostedTau_sv_x);
   fChain->SetBranchAddress("tau_sv_y", &tau_sv_y, &b_tau_sv_y);
   fChain->SetBranchAddress("boostedTau_sv_y", &boostedTau_sv_y, &b_boostedTau_sv_y);
   fChain->SetBranchAddress("tau_sv_z", &tau_sv_z, &b_tau_sv_z);
   fChain->SetBranchAddress("boostedTau_sv_z", &boostedTau_sv_z, &b_boostedTau_sv_z);
   fChain->SetBranchAddress("tau_flightLength_x", &tau_flightLength_x, &b_tau_flightLength_x);
   fChain->SetBranchAddress("boostedTau_flightLength_x", &boostedTau_flightLength_x, &b_boostedTau_flightLength_x);
   fChain->SetBranchAddress("tau_flightLength_y", &tau_flightLength_y, &b_tau_flightLength_y);
   fChain->SetBranchAddress("boostedTau_flightLength_y", &boostedTau_flightLength_y, &b_boostedTau_flightLength_y);
   fChain->SetBranchAddress("tau_flightLength_z", &tau_flightLength_z, &b_tau_flightLength_z);
   fChain->SetBranchAddress("boostedTau_flightLength_z", &boostedTau_flightLength_z, &b_boostedTau_flightLength_z);
   fChain->SetBranchAddress("tau_flightLength_sig", &tau_flightLength_sig, &b_tau_flightLength_sig);
   fChain->SetBranchAddress("boostedTau_flightLength_sig", &boostedTau_flightLength_sig, &b_boostedTau_flightLength_sig);
   fChain->SetBranchAddress("tau_pt_weighted_deta_strip", &tau_pt_weighted_deta_strip, &b_tau_pt_weighted_deta_strip);
   fChain->SetBranchAddress("boostedTau_pt_weighted_deta_strip", &boostedTau_pt_weighted_deta_strip, &b_boostedTau_pt_weighted_deta_strip);
   fChain->SetBranchAddress("tau_pt_weighted_dphi_strip", &tau_pt_weighted_dphi_strip, &b_tau_pt_weighted_dphi_strip);
   fChain->SetBranchAddress("boostedTau_pt_weighted_dphi_strip", &boostedTau_pt_weighted_dphi_strip, &b_boostedTau_pt_weighted_dphi_strip);
   fChain->SetBranchAddress("tau_pt_weighted_dr_signal", &tau_pt_weighted_dr_signal, &b_tau_pt_weighted_dr_signal);
   fChain->SetBranchAddress("boostedTau_pt_weighted_dr_signal", &boostedTau_pt_weighted_dr_signal, &b_boostedTau_pt_weighted_dr_signal);
   fChain->SetBranchAddress("tau_pt_weighted_dr_iso", &tau_pt_weighted_dr_iso, &b_tau_pt_weighted_dr_iso);
   fChain->SetBranchAddress("boostedTau_pt_weighted_dr_iso", &boostedTau_pt_weighted_dr_iso, &b_boostedTau_pt_weighted_dr_iso);
   fChain->SetBranchAddress("tau_leadingTrackNormChi2", &tau_leadingTrackNormChi2, &b_tau_leadingTrackNormChi2);
   fChain->SetBranchAddress("boostedTau_leadingTrackNormChi2", &boostedTau_leadingTrackNormChi2, &b_boostedTau_leadingTrackNormChi2);
   fChain->SetBranchAddress("tau_e_ratio", &tau_e_ratio, &b_tau_e_ratio);
   fChain->SetBranchAddress("boostedTau_e_ratio", &boostedTau_e_ratio, &b_boostedTau_e_ratio);
   fChain->SetBranchAddress("tau_gj_angle_diff", &tau_gj_angle_diff, &b_tau_gj_angle_diff);
   fChain->SetBranchAddress("boostedTau_gj_angle_diff", &boostedTau_gj_angle_diff, &b_boostedTau_gj_angle_diff);
   fChain->SetBranchAddress("tau_n_photons", &tau_n_photons, &b_tau_n_photons);
   fChain->SetBranchAddress("boostedTau_n_photons", &boostedTau_n_photons, &b_boostedTau_n_photons);
   fChain->SetBranchAddress("tau_emFraction", &tau_emFraction, &b_tau_emFraction);
   fChain->SetBranchAddress("boostedTau_emFraction", &boostedTau_emFraction, &b_boostedTau_emFraction);
   fChain->SetBranchAddress("tau_inside_ecal_crack", &tau_inside_ecal_crack, &b_tau_inside_ecal_crack);
   fChain->SetBranchAddress("boostedTau_inside_ecal_crack", &boostedTau_inside_ecal_crack, &b_boostedTau_inside_ecal_crack);
   fChain->SetBranchAddress("tau_leadChargedCand_etaAtEcalEntrance", &tau_leadChargedCand_etaAtEcalEntrance, &b_tau_leadChargedCand_etaAtEcalEntrance);
   fChain->SetBranchAddress("boostedTau_leadChargedCand_etaAtEcalEntrance", &boostedTau_leadChargedCand_etaAtEcalEntrance, &b_boostedTau_leadChargedCand_etaAtEcalEntrance);
   fChain->SetBranchAddress("pfCand_index", &pfCand_index, &b_pfCand_index);
   fChain->SetBranchAddress("lostTrack_index", &lostTrack_index, &b_lostTrack_index);
   fChain->SetBranchAddress("pfCand_tauSignal", &pfCand_tauSignal, &b_pfCand_tauSignal);
   fChain->SetBranchAddress("lostTrack_tauSignal", &lostTrack_tauSignal, &b_lostTrack_tauSignal);
   fChain->SetBranchAddress("pfCand_tauLeadChargedHadrCand", &pfCand_tauLeadChargedHadrCand, &b_pfCand_tauLeadChargedHadrCand);
   fChain->SetBranchAddress("lostTrack_tauLeadChargedHadrCand", &lostTrack_tauLeadChargedHadrCand, &b_lostTrack_tauLeadChargedHadrCand);
   fChain->SetBranchAddress("pfCand_tauIso", &pfCand_tauIso, &b_pfCand_tauIso);
   fChain->SetBranchAddress("lostTrack_tauIso", &lostTrack_tauIso, &b_lostTrack_tauIso);
   fChain->SetBranchAddress("pfCand_boostedTauSignal", &pfCand_boostedTauSignal, &b_pfCand_boostedTauSignal);
   fChain->SetBranchAddress("lostTrack_boostedTauSignal", &lostTrack_boostedTauSignal, &b_lostTrack_boostedTauSignal);
   fChain->SetBranchAddress("pfCand_boostedTauLeadChargedHadrCand", &pfCand_boostedTauLeadChargedHadrCand, &b_pfCand_boostedTauLeadChargedHadrCand);
   fChain->SetBranchAddress("lostTrack_boostedTauLeadChargedHadrCand", &lostTrack_boostedTauLeadChargedHadrCand, &b_lostTrack_boostedTauLeadChargedHadrCand);
   fChain->SetBranchAddress("pfCand_boostedTauIso", &pfCand_boostedTauIso, &b_pfCand_boostedTauIso);
   fChain->SetBranchAddress("lostTrack_boostedTauIso", &lostTrack_boostedTauIso, &b_lostTrack_boostedTauIso);
   fChain->SetBranchAddress("pfCand_jetDaughter", &pfCand_jetDaughter, &b_pfCand_jetDaughter);
   fChain->SetBranchAddress("lostTrack_jetDaughter", &lostTrack_jetDaughter, &b_lostTrack_jetDaughter);
   fChain->SetBranchAddress("pfCand_fatJetDaughter", &pfCand_fatJetDaughter, &b_pfCand_fatJetDaughter);
   fChain->SetBranchAddress("lostTrack_fatJetDaughter", &lostTrack_fatJetDaughter, &b_lostTrack_fatJetDaughter);
   fChain->SetBranchAddress("pfCand_subJetDaughter", &pfCand_subJetDaughter, &b_pfCand_subJetDaughter);
   fChain->SetBranchAddress("lostTrack_subJetDaughter", &lostTrack_subJetDaughter, &b_lostTrack_subJetDaughter);
   fChain->SetBranchAddress("pfCand_pt", &pfCand_pt, &b_pfCand_pt);
   fChain->SetBranchAddress("lostTrack_pt", &lostTrack_pt, &b_lostTrack_pt);
   fChain->SetBranchAddress("pfCand_eta", &pfCand_eta, &b_pfCand_eta);
   fChain->SetBranchAddress("lostTrack_eta", &lostTrack_eta, &b_lostTrack_eta);
   fChain->SetBranchAddress("pfCand_phi", &pfCand_phi, &b_pfCand_phi);
   fChain->SetBranchAddress("lostTrack_phi", &lostTrack_phi, &b_lostTrack_phi);
   fChain->SetBranchAddress("pfCand_mass", &pfCand_mass, &b_pfCand_mass);
   fChain->SetBranchAddress("lostTrack_mass", &lostTrack_mass, &b_lostTrack_mass);
   fChain->SetBranchAddress("pfCand_pvAssociationQuality", &pfCand_pvAssociationQuality, &b_pfCand_pvAssociationQuality);
   fChain->SetBranchAddress("lostTrack_pvAssociationQuality", &lostTrack_pvAssociationQuality, &b_lostTrack_pvAssociationQuality);
   fChain->SetBranchAddress("pfCand_fromPV", &pfCand_fromPV, &b_pfCand_fromPV);
   fChain->SetBranchAddress("lostTrack_fromPV", &lostTrack_fromPV, &b_lostTrack_fromPV);
   fChain->SetBranchAddress("pfCand_puppiWeight", &pfCand_puppiWeight, &b_pfCand_puppiWeight);
   fChain->SetBranchAddress("lostTrack_puppiWeight", &lostTrack_puppiWeight, &b_lostTrack_puppiWeight);
   fChain->SetBranchAddress("pfCand_puppiWeightNoLep", &pfCand_puppiWeightNoLep, &b_pfCand_puppiWeightNoLep);
   fChain->SetBranchAddress("lostTrack_puppiWeightNoLep", &lostTrack_puppiWeightNoLep, &b_lostTrack_puppiWeightNoLep);
   fChain->SetBranchAddress("pfCand_particleType", &pfCand_particleType, &b_pfCand_particleType);
   fChain->SetBranchAddress("lostTrack_particleType", &lostTrack_particleType, &b_lostTrack_particleType);
   fChain->SetBranchAddress("pfCand_charge", &pfCand_charge, &b_pfCand_charge);
   fChain->SetBranchAddress("lostTrack_charge", &lostTrack_charge, &b_lostTrack_charge);
   fChain->SetBranchAddress("pfCand_lostInnerHits", &pfCand_lostInnerHits, &b_pfCand_lostInnerHits);
   fChain->SetBranchAddress("lostTrack_lostInnerHits", &lostTrack_lostInnerHits, &b_lostTrack_lostInnerHits);
   fChain->SetBranchAddress("pfCand_nHits", &pfCand_nHits, &b_pfCand_nHits);
   fChain->SetBranchAddress("lostTrack_nHits", &lostTrack_nHits, &b_lostTrack_nHits);
   fChain->SetBranchAddress("pfCand_nPixelHits", &pfCand_nPixelHits, &b_pfCand_nPixelHits);
   fChain->SetBranchAddress("lostTrack_nPixelHits", &lostTrack_nPixelHits, &b_lostTrack_nPixelHits);
   fChain->SetBranchAddress("pfCand_nPixelLayers", &pfCand_nPixelLayers, &b_pfCand_nPixelLayers);
   fChain->SetBranchAddress("lostTrack_nPixelLayers", &lostTrack_nPixelLayers, &b_lostTrack_nPixelLayers);
   fChain->SetBranchAddress("pfCand_nStripLayers", &pfCand_nStripLayers, &b_pfCand_nStripLayers);
   fChain->SetBranchAddress("lostTrack_nStripLayers", &lostTrack_nStripLayers, &b_lostTrack_nStripLayers);
   fChain->SetBranchAddress("pfCand_vertex_x", &pfCand_vertex_x, &b_pfCand_vertex_x);
   fChain->SetBranchAddress("lostTrack_vertex_x", &lostTrack_vertex_x, &b_lostTrack_vertex_x);
   fChain->SetBranchAddress("pfCand_vertex_y", &pfCand_vertex_y, &b_pfCand_vertex_y);
   fChain->SetBranchAddress("lostTrack_vertex_y", &lostTrack_vertex_y, &b_lostTrack_vertex_y);
   fChain->SetBranchAddress("pfCand_vertex_z", &pfCand_vertex_z, &b_pfCand_vertex_z);
   fChain->SetBranchAddress("lostTrack_vertex_z", &lostTrack_vertex_z, &b_lostTrack_vertex_z);
   fChain->SetBranchAddress("pfCand_vertex_t", &pfCand_vertex_t, &b_pfCand_vertex_t);
   fChain->SetBranchAddress("lostTrack_vertex_t", &lostTrack_vertex_t, &b_lostTrack_vertex_t);
   fChain->SetBranchAddress("pfCand_time", &pfCand_time, &b_pfCand_time);
   fChain->SetBranchAddress("lostTrack_time", &lostTrack_time, &b_lostTrack_time);
   fChain->SetBranchAddress("pfCand_timeError", &pfCand_timeError, &b_pfCand_timeError);
   fChain->SetBranchAddress("lostTrack_timeError", &lostTrack_timeError, &b_lostTrack_timeError);
   fChain->SetBranchAddress("pfCand_hasTrackDetails", &pfCand_hasTrackDetails, &b_pfCand_hasTrackDetails);
   fChain->SetBranchAddress("lostTrack_hasTrackDetails", &lostTrack_hasTrackDetails, &b_lostTrack_hasTrackDetails);
   fChain->SetBranchAddress("pfCand_dxy", &pfCand_dxy, &b_pfCand_dxy);
   fChain->SetBranchAddress("lostTrack_dxy", &lostTrack_dxy, &b_lostTrack_dxy);
   fChain->SetBranchAddress("pfCand_dxy_error", &pfCand_dxy_error, &b_pfCand_dxy_error);
   fChain->SetBranchAddress("lostTrack_dxy_error", &lostTrack_dxy_error, &b_lostTrack_dxy_error);
   fChain->SetBranchAddress("pfCand_dz", &pfCand_dz, &b_pfCand_dz);
   fChain->SetBranchAddress("lostTrack_dz", &lostTrack_dz, &b_lostTrack_dz);
   fChain->SetBranchAddress("pfCand_dz_error", &pfCand_dz_error, &b_pfCand_dz_error);
   fChain->SetBranchAddress("lostTrack_dz_error", &lostTrack_dz_error, &b_lostTrack_dz_error);
   fChain->SetBranchAddress("pfCand_track_pt", &pfCand_track_pt, &b_pfCand_track_pt);
   fChain->SetBranchAddress("lostTrack_track_pt", &lostTrack_track_pt, &b_lostTrack_track_pt);
   fChain->SetBranchAddress("pfCand_track_eta", &pfCand_track_eta, &b_pfCand_track_eta);
   fChain->SetBranchAddress("lostTrack_track_eta", &lostTrack_track_eta, &b_lostTrack_track_eta);
   fChain->SetBranchAddress("pfCand_track_phi", &pfCand_track_phi, &b_pfCand_track_phi);
   fChain->SetBranchAddress("lostTrack_track_phi", &lostTrack_track_phi, &b_lostTrack_track_phi);
   fChain->SetBranchAddress("pfCand_track_chi2", &pfCand_track_chi2, &b_pfCand_track_chi2);
   fChain->SetBranchAddress("lostTrack_track_chi2", &lostTrack_track_chi2, &b_lostTrack_track_chi2);
   fChain->SetBranchAddress("pfCand_track_ndof", &pfCand_track_ndof, &b_pfCand_track_ndof);
   fChain->SetBranchAddress("lostTrack_track_ndof", &lostTrack_track_ndof, &b_lostTrack_track_ndof);
   fChain->SetBranchAddress("pfCand_caloFraction", &pfCand_caloFraction, &b_pfCand_caloFraction);
   fChain->SetBranchAddress("lostTrack_caloFraction", &lostTrack_caloFraction, &b_lostTrack_caloFraction);
   fChain->SetBranchAddress("pfCand_hcalFraction", &pfCand_hcalFraction, &b_pfCand_hcalFraction);
   fChain->SetBranchAddress("lostTrack_hcalFraction", &lostTrack_hcalFraction, &b_lostTrack_hcalFraction);
   fChain->SetBranchAddress("pfCand_rawCaloFraction", &pfCand_rawCaloFraction, &b_pfCand_rawCaloFraction);
   fChain->SetBranchAddress("lostTrack_rawCaloFraction", &lostTrack_rawCaloFraction, &b_lostTrack_rawCaloFraction);
   fChain->SetBranchAddress("pfCand_rawHcalFraction", &pfCand_rawHcalFraction, &b_pfCand_rawHcalFraction);
   fChain->SetBranchAddress("lostTrack_rawHcalFraction", &lostTrack_rawHcalFraction, &b_lostTrack_rawHcalFraction);
   fChain->SetBranchAddress("ele_index", &ele_index, &b_ele_index);
   fChain->SetBranchAddress("ele_pt", &ele_pt, &b_ele_pt);
   fChain->SetBranchAddress("ele_eta", &ele_eta, &b_ele_eta);
   fChain->SetBranchAddress("ele_phi", &ele_phi, &b_ele_phi);
   fChain->SetBranchAddress("ele_mass", &ele_mass, &b_ele_mass);
   fChain->SetBranchAddress("ele_cc_ele_energy", &ele_cc_ele_energy, &b_ele_cc_ele_energy);
   fChain->SetBranchAddress("ele_cc_gamma_energy", &ele_cc_gamma_energy, &b_ele_cc_gamma_energy);
   fChain->SetBranchAddress("ele_cc_n_gamma", &ele_cc_n_gamma, &b_ele_cc_n_gamma);
   fChain->SetBranchAddress("ele_dxy", &ele_dxy, &b_ele_dxy);
   fChain->SetBranchAddress("ele_dxy_error", &ele_dxy_error, &b_ele_dxy_error);
   fChain->SetBranchAddress("ele_ip3d", &ele_ip3d, &b_ele_ip3d);
   fChain->SetBranchAddress("ele_trackMomentumAtVtx", &ele_trackMomentumAtVtx, &b_ele_trackMomentumAtVtx);
   fChain->SetBranchAddress("ele_trackMomentumAtCalo", &ele_trackMomentumAtCalo, &b_ele_trackMomentumAtCalo);
   fChain->SetBranchAddress("ele_trackMomentumOut", &ele_trackMomentumOut, &b_ele_trackMomentumOut);
   fChain->SetBranchAddress("ele_trackMomentumAtEleClus", &ele_trackMomentumAtEleClus, &b_ele_trackMomentumAtEleClus);
   fChain->SetBranchAddress("ele_trackMomentumAtVtxWithConstraint", &ele_trackMomentumAtVtxWithConstraint, &b_ele_trackMomentumAtVtxWithConstraint);
   fChain->SetBranchAddress("ele_ecalEnergy", &ele_ecalEnergy, &b_ele_ecalEnergy);
   fChain->SetBranchAddress("ele_ecalEnergy_error", &ele_ecalEnergy_error, &b_ele_ecalEnergy_error);
   fChain->SetBranchAddress("ele_eSuperClusterOverP", &ele_eSuperClusterOverP, &b_ele_eSuperClusterOverP);
   fChain->SetBranchAddress("ele_eSeedClusterOverP", &ele_eSeedClusterOverP, &b_ele_eSeedClusterOverP);
   fChain->SetBranchAddress("ele_eSeedClusterOverPout", &ele_eSeedClusterOverPout, &b_ele_eSeedClusterOverPout);
   fChain->SetBranchAddress("ele_eEleClusterOverPout", &ele_eEleClusterOverPout, &b_ele_eEleClusterOverPout);
   fChain->SetBranchAddress("ele_deltaEtaSuperClusterTrackAtVtx", &ele_deltaEtaSuperClusterTrackAtVtx, &b_ele_deltaEtaSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("ele_deltaEtaSeedClusterTrackAtCalo", &ele_deltaEtaSeedClusterTrackAtCalo, &b_ele_deltaEtaSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("ele_deltaEtaEleClusterTrackAtCalo", &ele_deltaEtaEleClusterTrackAtCalo, &b_ele_deltaEtaEleClusterTrackAtCalo);
   fChain->SetBranchAddress("ele_deltaEtaSeedClusterTrackAtVtx", &ele_deltaEtaSeedClusterTrackAtVtx, &b_ele_deltaEtaSeedClusterTrackAtVtx);
   fChain->SetBranchAddress("ele_deltaPhiEleClusterTrackAtCalo", &ele_deltaPhiEleClusterTrackAtCalo, &b_ele_deltaPhiEleClusterTrackAtCalo);
   fChain->SetBranchAddress("ele_deltaPhiSuperClusterTrackAtVtx", &ele_deltaPhiSuperClusterTrackAtVtx, &b_ele_deltaPhiSuperClusterTrackAtVtx);
   fChain->SetBranchAddress("ele_deltaPhiSeedClusterTrackAtCalo", &ele_deltaPhiSeedClusterTrackAtCalo, &b_ele_deltaPhiSeedClusterTrackAtCalo);
   fChain->SetBranchAddress("ele_mvaInput_earlyBrem", &ele_mvaInput_earlyBrem, &b_ele_mvaInput_earlyBrem);
   fChain->SetBranchAddress("ele_mvaInput_lateBrem", &ele_mvaInput_lateBrem, &b_ele_mvaInput_lateBrem);
   fChain->SetBranchAddress("ele_mvaInput_sigmaEtaEta", &ele_mvaInput_sigmaEtaEta, &b_ele_mvaInput_sigmaEtaEta);
   fChain->SetBranchAddress("ele_mvaInput_hadEnergy", &ele_mvaInput_hadEnergy, &b_ele_mvaInput_hadEnergy);
   fChain->SetBranchAddress("ele_mvaInput_deltaEta", &ele_mvaInput_deltaEta, &b_ele_mvaInput_deltaEta);
   fChain->SetBranchAddress("ele_gsfTrack_normalizedChi2", &ele_gsfTrack_normalizedChi2, &b_ele_gsfTrack_normalizedChi2);
   fChain->SetBranchAddress("ele_gsfTrack_numberOfValidHits", &ele_gsfTrack_numberOfValidHits, &b_ele_gsfTrack_numberOfValidHits);
   fChain->SetBranchAddress("ele_gsfTrack_pt", &ele_gsfTrack_pt, &b_ele_gsfTrack_pt);
   fChain->SetBranchAddress("ele_gsfTrack_pt_error", &ele_gsfTrack_pt_error, &b_ele_gsfTrack_pt_error);
   fChain->SetBranchAddress("ele_closestCtfTrack_normalizedChi2", &ele_closestCtfTrack_normalizedChi2, &b_ele_closestCtfTrack_normalizedChi2);
   fChain->SetBranchAddress("ele_closestCtfTrack_numberOfValidHits", &ele_closestCtfTrack_numberOfValidHits, &b_ele_closestCtfTrack_numberOfValidHits);
   fChain->SetBranchAddress("ele_sigmaEtaEta", &ele_sigmaEtaEta, &b_ele_sigmaEtaEta);
   fChain->SetBranchAddress("ele_sigmaIetaIeta", &ele_sigmaIetaIeta, &b_ele_sigmaIetaIeta);
   fChain->SetBranchAddress("ele_sigmaIphiIphi", &ele_sigmaIphiIphi, &b_ele_sigmaIphiIphi);
   fChain->SetBranchAddress("ele_sigmaIetaIphi", &ele_sigmaIetaIphi, &b_ele_sigmaIetaIphi);
   fChain->SetBranchAddress("ele_e1x5", &ele_e1x5, &b_ele_e1x5);
   fChain->SetBranchAddress("ele_e2x5Max", &ele_e2x5Max, &b_ele_e2x5Max);
   fChain->SetBranchAddress("ele_e5x5", &ele_e5x5, &b_ele_e5x5);
   fChain->SetBranchAddress("ele_r9", &ele_r9, &b_ele_r9);
   fChain->SetBranchAddress("ele_hcalDepth1OverEcal", &ele_hcalDepth1OverEcal, &b_ele_hcalDepth1OverEcal);
   fChain->SetBranchAddress("ele_hcalDepth2OverEcal", &ele_hcalDepth2OverEcal, &b_ele_hcalDepth2OverEcal);
   fChain->SetBranchAddress("ele_hcalDepth1OverEcalBc", &ele_hcalDepth1OverEcalBc, &b_ele_hcalDepth1OverEcalBc);
   fChain->SetBranchAddress("ele_hcalDepth2OverEcalBc", &ele_hcalDepth2OverEcalBc, &b_ele_hcalDepth2OverEcalBc);
   fChain->SetBranchAddress("ele_eLeft", &ele_eLeft, &b_ele_eLeft);
   fChain->SetBranchAddress("ele_eRight", &ele_eRight, &b_ele_eRight);
   fChain->SetBranchAddress("ele_eBottom", &ele_eBottom, &b_ele_eBottom);
   fChain->SetBranchAddress("ele_eTop", &ele_eTop, &b_ele_eTop);
   fChain->SetBranchAddress("ele_full5x5_sigmaEtaEta", &ele_full5x5_sigmaEtaEta, &b_ele_full5x5_sigmaEtaEta);
   fChain->SetBranchAddress("ele_full5x5_sigmaIetaIeta", &ele_full5x5_sigmaIetaIeta, &b_ele_full5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("ele_full5x5_sigmaIphiIphi", &ele_full5x5_sigmaIphiIphi, &b_ele_full5x5_sigmaIphiIphi);
   fChain->SetBranchAddress("ele_full5x5_sigmaIetaIphi", &ele_full5x5_sigmaIetaIphi, &b_ele_full5x5_sigmaIetaIphi);
   fChain->SetBranchAddress("ele_full5x5_e1x5", &ele_full5x5_e1x5, &b_ele_full5x5_e1x5);
   fChain->SetBranchAddress("ele_full5x5_e2x5Max", &ele_full5x5_e2x5Max, &b_ele_full5x5_e2x5Max);
   fChain->SetBranchAddress("ele_full5x5_e5x5", &ele_full5x5_e5x5, &b_ele_full5x5_e5x5);
   fChain->SetBranchAddress("ele_full5x5_r9", &ele_full5x5_r9, &b_ele_full5x5_r9);
   fChain->SetBranchAddress("ele_full5x5_hcalDepth1OverEcal", &ele_full5x5_hcalDepth1OverEcal, &b_ele_full5x5_hcalDepth1OverEcal);
   fChain->SetBranchAddress("ele_full5x5_hcalDepth2OverEcal", &ele_full5x5_hcalDepth2OverEcal, &b_ele_full5x5_hcalDepth2OverEcal);
   fChain->SetBranchAddress("ele_full5x5_hcalDepth1OverEcalBc", &ele_full5x5_hcalDepth1OverEcalBc, &b_ele_full5x5_hcalDepth1OverEcalBc);
   fChain->SetBranchAddress("ele_full5x5_hcalDepth2OverEcalBc", &ele_full5x5_hcalDepth2OverEcalBc, &b_ele_full5x5_hcalDepth2OverEcalBc);
   fChain->SetBranchAddress("ele_full5x5_eLeft", &ele_full5x5_eLeft, &b_ele_full5x5_eLeft);
   fChain->SetBranchAddress("ele_full5x5_eRight", &ele_full5x5_eRight, &b_ele_full5x5_eRight);
   fChain->SetBranchAddress("ele_full5x5_eBottom", &ele_full5x5_eBottom, &b_ele_full5x5_eBottom);
   fChain->SetBranchAddress("ele_full5x5_eTop", &ele_full5x5_eTop, &b_ele_full5x5_eTop);
   fChain->SetBranchAddress("ele_full5x5_e2x5Left", &ele_full5x5_e2x5Left, &b_ele_full5x5_e2x5Left);
   fChain->SetBranchAddress("ele_full5x5_e2x5Right", &ele_full5x5_e2x5Right, &b_ele_full5x5_e2x5Right);
   fChain->SetBranchAddress("ele_full5x5_e2x5Bottom", &ele_full5x5_e2x5Bottom, &b_ele_full5x5_e2x5Bottom);
   fChain->SetBranchAddress("ele_full5x5_e2x5Top", &ele_full5x5_e2x5Top, &b_ele_full5x5_e2x5Top);
   fChain->SetBranchAddress("ele_hgcal_sigmaUU", &ele_hgcal_sigmaUU, &b_ele_hgcal_sigmaUU);
   fChain->SetBranchAddress("ele_hgcal_sigmaVV", &ele_hgcal_sigmaVV, &b_ele_hgcal_sigmaVV);
   fChain->SetBranchAddress("ele_hgcal_sigmaEE", &ele_hgcal_sigmaEE, &b_ele_hgcal_sigmaEE);
   fChain->SetBranchAddress("ele_hgcal_sigmaPP", &ele_hgcal_sigmaPP, &b_ele_hgcal_sigmaPP);
   fChain->SetBranchAddress("ele_hgcal_nLayers", &ele_hgcal_nLayers, &b_ele_hgcal_nLayers);
   fChain->SetBranchAddress("ele_hgcal_firstLayer", &ele_hgcal_firstLayer, &b_ele_hgcal_firstLayer);
   fChain->SetBranchAddress("ele_hgcal_lastLayer", &ele_hgcal_lastLayer, &b_ele_hgcal_lastLayer);
   fChain->SetBranchAddress("ele_hgcal_layerEfrac10", &ele_hgcal_layerEfrac10, &b_ele_hgcal_layerEfrac10);
   fChain->SetBranchAddress("ele_hgcal_layerEfrac90", &ele_hgcal_layerEfrac90, &b_ele_hgcal_layerEfrac90);
   fChain->SetBranchAddress("ele_hgcal_e4oEtot", &ele_hgcal_e4oEtot, &b_ele_hgcal_e4oEtot);
   fChain->SetBranchAddress("ele_hgcal_ecEnergy", &ele_hgcal_ecEnergy, &b_ele_hgcal_ecEnergy);
   fChain->SetBranchAddress("ele_hgcal_ecEnergyEE", &ele_hgcal_ecEnergyEE, &b_ele_hgcal_ecEnergyEE);
   fChain->SetBranchAddress("ele_hgcal_ecEnergyFH", &ele_hgcal_ecEnergyFH, &b_ele_hgcal_ecEnergyFH);
   fChain->SetBranchAddress("ele_hgcal_ecEnergyBH", &ele_hgcal_ecEnergyBH, &b_ele_hgcal_ecEnergyBH);
   fChain->SetBranchAddress("ele_hgcal_ecEt", &ele_hgcal_ecEt, &b_ele_hgcal_ecEt);
   fChain->SetBranchAddress("ele_hgcal_ecOrigEnergy", &ele_hgcal_ecOrigEnergy, &b_ele_hgcal_ecOrigEnergy);
   fChain->SetBranchAddress("ele_hgcal_ecOrigEt", &ele_hgcal_ecOrigEt, &b_ele_hgcal_ecOrigEt);
   fChain->SetBranchAddress("ele_hgcal_caloIsoRing0", &ele_hgcal_caloIsoRing0, &b_ele_hgcal_caloIsoRing0);
   fChain->SetBranchAddress("ele_hgcal_caloIsoRing1", &ele_hgcal_caloIsoRing1, &b_ele_hgcal_caloIsoRing1);
   fChain->SetBranchAddress("ele_hgcal_caloIsoRing2", &ele_hgcal_caloIsoRing2, &b_ele_hgcal_caloIsoRing2);
   fChain->SetBranchAddress("ele_hgcal_caloIsoRing3", &ele_hgcal_caloIsoRing3, &b_ele_hgcal_caloIsoRing3);
   fChain->SetBranchAddress("ele_hgcal_caloIsoRing4", &ele_hgcal_caloIsoRing4, &b_ele_hgcal_caloIsoRing4);
   fChain->SetBranchAddress("ele_hgcal_depthCompatibility", &ele_hgcal_depthCompatibility, &b_ele_hgcal_depthCompatibility);
   fChain->SetBranchAddress("ele_hgcal_expectedDepth", &ele_hgcal_expectedDepth, &b_ele_hgcal_expectedDepth);
   fChain->SetBranchAddress("ele_hgcal_expectedSigma", &ele_hgcal_expectedSigma, &b_ele_hgcal_expectedSigma);
   fChain->SetBranchAddress("ele_hgcal_measuredDepth", &ele_hgcal_measuredDepth, &b_ele_hgcal_measuredDepth);
   fChain->SetBranchAddress("ele_hgcal_pcaAxisX", &ele_hgcal_pcaAxisX, &b_ele_hgcal_pcaAxisX);
   fChain->SetBranchAddress("ele_hgcal_pcaAxisY", &ele_hgcal_pcaAxisY, &b_ele_hgcal_pcaAxisY);
   fChain->SetBranchAddress("ele_hgcal_pcaAxisZ", &ele_hgcal_pcaAxisZ, &b_ele_hgcal_pcaAxisZ);
   fChain->SetBranchAddress("ele_hgcal_pcaPositionX", &ele_hgcal_pcaPositionX, &b_ele_hgcal_pcaPositionX);
   fChain->SetBranchAddress("ele_hgcal_pcaPositionY", &ele_hgcal_pcaPositionY, &b_ele_hgcal_pcaPositionY);
   fChain->SetBranchAddress("ele_hgcal_pcaPositionZ", &ele_hgcal_pcaPositionZ, &b_ele_hgcal_pcaPositionZ);
   fChain->SetBranchAddress("ele_hgcal_pcaEig1", &ele_hgcal_pcaEig1, &b_ele_hgcal_pcaEig1);
   fChain->SetBranchAddress("ele_hgcal_pcaEig2", &ele_hgcal_pcaEig2, &b_ele_hgcal_pcaEig2);
   fChain->SetBranchAddress("ele_hgcal_pcaEig3", &ele_hgcal_pcaEig3, &b_ele_hgcal_pcaEig3);
   fChain->SetBranchAddress("ele_hgcal_pcaSig1", &ele_hgcal_pcaSig1, &b_ele_hgcal_pcaSig1);
   fChain->SetBranchAddress("ele_hgcal_pcaSig2", &ele_hgcal_pcaSig2, &b_ele_hgcal_pcaSig2);
   fChain->SetBranchAddress("ele_hgcal_pcaSig3", &ele_hgcal_pcaSig3, &b_ele_hgcal_pcaSig3);
   fChain->SetBranchAddress("muon_index", &muon_index, &b_muon_index);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_mass", &muon_mass, &b_muon_mass);
   fChain->SetBranchAddress("muon_dxy", &muon_dxy, &b_muon_dxy);
   fChain->SetBranchAddress("muon_dxy_error", &muon_dxy_error, &b_muon_dxy_error);
   fChain->SetBranchAddress("muon_normalizedChi2", &muon_normalizedChi2, &b_muon_normalizedChi2);
   fChain->SetBranchAddress("muon_numberOfValidHits", &muon_numberOfValidHits, &b_muon_numberOfValidHits);
   fChain->SetBranchAddress("muon_segmentCompatibility", &muon_segmentCompatibility, &b_muon_segmentCompatibility);
   fChain->SetBranchAddress("muon_caloCompatibility", &muon_caloCompatibility, &b_muon_caloCompatibility);
   fChain->SetBranchAddress("muon_pfEcalEnergy", &muon_pfEcalEnergy, &b_muon_pfEcalEnergy);
   fChain->SetBranchAddress("muon_type", &muon_type, &b_muon_type);
   fChain->SetBranchAddress("muon_n_matches_DT_1", &muon_n_matches_DT_1, &b_muon_n_matches_DT_1);
   fChain->SetBranchAddress("muon_n_matches_DT_2", &muon_n_matches_DT_2, &b_muon_n_matches_DT_2);
   fChain->SetBranchAddress("muon_n_matches_DT_3", &muon_n_matches_DT_3, &b_muon_n_matches_DT_3);
   fChain->SetBranchAddress("muon_n_matches_DT_4", &muon_n_matches_DT_4, &b_muon_n_matches_DT_4);
   fChain->SetBranchAddress("muon_n_matches_CSC_1", &muon_n_matches_CSC_1, &b_muon_n_matches_CSC_1);
   fChain->SetBranchAddress("muon_n_matches_CSC_2", &muon_n_matches_CSC_2, &b_muon_n_matches_CSC_2);
   fChain->SetBranchAddress("muon_n_matches_CSC_3", &muon_n_matches_CSC_3, &b_muon_n_matches_CSC_3);
   fChain->SetBranchAddress("muon_n_matches_CSC_4", &muon_n_matches_CSC_4, &b_muon_n_matches_CSC_4);
   fChain->SetBranchAddress("muon_n_matches_RPC_1", &muon_n_matches_RPC_1, &b_muon_n_matches_RPC_1);
   fChain->SetBranchAddress("muon_n_matches_RPC_2", &muon_n_matches_RPC_2, &b_muon_n_matches_RPC_2);
   fChain->SetBranchAddress("muon_n_matches_RPC_3", &muon_n_matches_RPC_3, &b_muon_n_matches_RPC_3);
   fChain->SetBranchAddress("muon_n_matches_RPC_4", &muon_n_matches_RPC_4, &b_muon_n_matches_RPC_4);
   fChain->SetBranchAddress("muon_n_matches_GEM_1", &muon_n_matches_GEM_1, &b_muon_n_matches_GEM_1);
   fChain->SetBranchAddress("muon_n_matches_GEM_2", &muon_n_matches_GEM_2, &b_muon_n_matches_GEM_2);
   fChain->SetBranchAddress("muon_n_matches_GEM_3", &muon_n_matches_GEM_3, &b_muon_n_matches_GEM_3);
   fChain->SetBranchAddress("muon_n_matches_GEM_4", &muon_n_matches_GEM_4, &b_muon_n_matches_GEM_4);
   fChain->SetBranchAddress("muon_n_matches_ME0_1", &muon_n_matches_ME0_1, &b_muon_n_matches_ME0_1);
   fChain->SetBranchAddress("muon_n_matches_ME0_2", &muon_n_matches_ME0_2, &b_muon_n_matches_ME0_2);
   fChain->SetBranchAddress("muon_n_matches_ME0_3", &muon_n_matches_ME0_3, &b_muon_n_matches_ME0_3);
   fChain->SetBranchAddress("muon_n_matches_ME0_4", &muon_n_matches_ME0_4, &b_muon_n_matches_ME0_4);
   fChain->SetBranchAddress("muon_n_hits_DT_1", &muon_n_hits_DT_1, &b_muon_n_hits_DT_1);
   fChain->SetBranchAddress("muon_n_hits_DT_2", &muon_n_hits_DT_2, &b_muon_n_hits_DT_2);
   fChain->SetBranchAddress("muon_n_hits_DT_3", &muon_n_hits_DT_3, &b_muon_n_hits_DT_3);
   fChain->SetBranchAddress("muon_n_hits_DT_4", &muon_n_hits_DT_4, &b_muon_n_hits_DT_4);
   fChain->SetBranchAddress("muon_n_hits_CSC_1", &muon_n_hits_CSC_1, &b_muon_n_hits_CSC_1);
   fChain->SetBranchAddress("muon_n_hits_CSC_2", &muon_n_hits_CSC_2, &b_muon_n_hits_CSC_2);
   fChain->SetBranchAddress("muon_n_hits_CSC_3", &muon_n_hits_CSC_3, &b_muon_n_hits_CSC_3);
   fChain->SetBranchAddress("muon_n_hits_CSC_4", &muon_n_hits_CSC_4, &b_muon_n_hits_CSC_4);
   fChain->SetBranchAddress("muon_n_hits_RPC_1", &muon_n_hits_RPC_1, &b_muon_n_hits_RPC_1);
   fChain->SetBranchAddress("muon_n_hits_RPC_2", &muon_n_hits_RPC_2, &b_muon_n_hits_RPC_2);
   fChain->SetBranchAddress("muon_n_hits_RPC_3", &muon_n_hits_RPC_3, &b_muon_n_hits_RPC_3);
   fChain->SetBranchAddress("muon_n_hits_RPC_4", &muon_n_hits_RPC_4, &b_muon_n_hits_RPC_4);
   fChain->SetBranchAddress("muon_n_hits_GEM_1", &muon_n_hits_GEM_1, &b_muon_n_hits_GEM_1);
   fChain->SetBranchAddress("muon_n_hits_GEM_2", &muon_n_hits_GEM_2, &b_muon_n_hits_GEM_2);
   fChain->SetBranchAddress("muon_n_hits_GEM_3", &muon_n_hits_GEM_3, &b_muon_n_hits_GEM_3);
   fChain->SetBranchAddress("muon_n_hits_GEM_4", &muon_n_hits_GEM_4, &b_muon_n_hits_GEM_4);
   fChain->SetBranchAddress("muon_n_hits_ME0_1", &muon_n_hits_ME0_1, &b_muon_n_hits_ME0_1);
   fChain->SetBranchAddress("muon_n_hits_ME0_2", &muon_n_hits_ME0_2, &b_muon_n_hits_ME0_2);
   fChain->SetBranchAddress("muon_n_hits_ME0_3", &muon_n_hits_ME0_3, &b_muon_n_hits_ME0_3);
   fChain->SetBranchAddress("muon_n_hits_ME0_4", &muon_n_hits_ME0_4, &b_muon_n_hits_ME0_4);
   fChain->SetBranchAddress("isoTrack_index", &isoTrack_index, &b_isoTrack_index);
   fChain->SetBranchAddress("isoTrack_pt", &isoTrack_pt, &b_isoTrack_pt);
   fChain->SetBranchAddress("isoTrack_eta", &isoTrack_eta, &b_isoTrack_eta);
   fChain->SetBranchAddress("isoTrack_phi", &isoTrack_phi, &b_isoTrack_phi);
   fChain->SetBranchAddress("isoTrack_fromPV", &isoTrack_fromPV, &b_isoTrack_fromPV);
   fChain->SetBranchAddress("isoTrack_charge", &isoTrack_charge, &b_isoTrack_charge);
   fChain->SetBranchAddress("isoTrack_dxy", &isoTrack_dxy, &b_isoTrack_dxy);
   fChain->SetBranchAddress("isoTrack_dxy_error", &isoTrack_dxy_error, &b_isoTrack_dxy_error);
   fChain->SetBranchAddress("isoTrack_dz", &isoTrack_dz, &b_isoTrack_dz);
   fChain->SetBranchAddress("isoTrack_dz_error", &isoTrack_dz_error, &b_isoTrack_dz_error);
   fChain->SetBranchAddress("isoTrack_isHighPurityTrack", &isoTrack_isHighPurityTrack, &b_isoTrack_isHighPurityTrack);
   fChain->SetBranchAddress("isoTrack_isTightTrack", &isoTrack_isTightTrack, &b_isoTrack_isTightTrack);
   fChain->SetBranchAddress("isoTrack_isLooseTrack", &isoTrack_isLooseTrack, &b_isoTrack_isLooseTrack);
   fChain->SetBranchAddress("isoTrack_dEdxStrip", &isoTrack_dEdxStrip, &b_isoTrack_dEdxStrip);
   fChain->SetBranchAddress("isoTrack_dEdxPixel", &isoTrack_dEdxPixel, &b_isoTrack_dEdxPixel);
   fChain->SetBranchAddress("isoTrack_deltaEta", &isoTrack_deltaEta, &b_isoTrack_deltaEta);
   fChain->SetBranchAddress("isoTrack_deltaPhi", &isoTrack_deltaPhi, &b_isoTrack_deltaPhi);
   fChain->SetBranchAddress("isoTrack_n_ValidHits", &isoTrack_n_ValidHits, &b_isoTrack_n_ValidHits);
   fChain->SetBranchAddress("isoTrack_n_BadHits", &isoTrack_n_BadHits, &b_isoTrack_n_BadHits);
   fChain->SetBranchAddress("isoTrack_n_TimingHits", &isoTrack_n_TimingHits, &b_isoTrack_n_TimingHits);
   fChain->SetBranchAddress("isoTrack_n_ValidTimingHits", &isoTrack_n_ValidTimingHits, &b_isoTrack_n_ValidTimingHits);
   fChain->SetBranchAddress("isoTrack_n_LostTimingHits", &isoTrack_n_LostTimingHits, &b_isoTrack_n_LostTimingHits);
   fChain->SetBranchAddress("isoTrack_n_MuonHits", &isoTrack_n_MuonHits, &b_isoTrack_n_MuonHits);
   fChain->SetBranchAddress("isoTrack_n_ValidMuonHits", &isoTrack_n_ValidMuonHits, &b_isoTrack_n_ValidMuonHits);
   fChain->SetBranchAddress("isoTrack_n_LostMuonHits", &isoTrack_n_LostMuonHits, &b_isoTrack_n_LostMuonHits);
   fChain->SetBranchAddress("isoTrack_n_BadMuonHits", &isoTrack_n_BadMuonHits, &b_isoTrack_n_BadMuonHits);
   fChain->SetBranchAddress("isoTrack_n_ValidMuonDTHits", &isoTrack_n_ValidMuonDTHits, &b_isoTrack_n_ValidMuonDTHits);
   fChain->SetBranchAddress("isoTrack_n_LostMuonDTHits", &isoTrack_n_LostMuonDTHits, &b_isoTrack_n_LostMuonDTHits);
   fChain->SetBranchAddress("isoTrack_n_BadMuonDTHits", &isoTrack_n_BadMuonDTHits, &b_isoTrack_n_BadMuonDTHits);
   fChain->SetBranchAddress("isoTrack_n_ValidMuonCSCHits", &isoTrack_n_ValidMuonCSCHits, &b_isoTrack_n_ValidMuonCSCHits);
   fChain->SetBranchAddress("isoTrack_n_LostMuonCSCHits", &isoTrack_n_LostMuonCSCHits, &b_isoTrack_n_LostMuonCSCHits);
   fChain->SetBranchAddress("isoTrack_n_BadMuonCSCHits", &isoTrack_n_BadMuonCSCHits, &b_isoTrack_n_BadMuonCSCHits);
   fChain->SetBranchAddress("isoTrack_n_ValidMuonRPCHits", &isoTrack_n_ValidMuonRPCHits, &b_isoTrack_n_ValidMuonRPCHits);
   fChain->SetBranchAddress("isoTrack_n_LostMuonRPCHits", &isoTrack_n_LostMuonRPCHits, &b_isoTrack_n_LostMuonRPCHits);
   fChain->SetBranchAddress("isoTrack_n_BadMuonRPCHits", &isoTrack_n_BadMuonRPCHits, &b_isoTrack_n_BadMuonRPCHits);
   fChain->SetBranchAddress("isoTrack_n_ValidMuonGEMHits", &isoTrack_n_ValidMuonGEMHits, &b_isoTrack_n_ValidMuonGEMHits);
   fChain->SetBranchAddress("isoTrack_n_LostMuonGEMHits", &isoTrack_n_LostMuonGEMHits, &b_isoTrack_n_LostMuonGEMHits);
   fChain->SetBranchAddress("isoTrack_n_BadMuonGEMHits", &isoTrack_n_BadMuonGEMHits, &b_isoTrack_n_BadMuonGEMHits);
   fChain->SetBranchAddress("isoTrack_n_ValidMuonME0Hits", &isoTrack_n_ValidMuonME0Hits, &b_isoTrack_n_ValidMuonME0Hits);
   fChain->SetBranchAddress("isoTrack_n_LostMuonME0Hits", &isoTrack_n_LostMuonME0Hits, &b_isoTrack_n_LostMuonME0Hits);
   fChain->SetBranchAddress("isoTrack_n_BadMuonME0Hits", &isoTrack_n_BadMuonME0Hits, &b_isoTrack_n_BadMuonME0Hits);
   fChain->SetBranchAddress("isoTrack_n_InactiveHits", &isoTrack_n_InactiveHits, &b_isoTrack_n_InactiveHits);
   fChain->SetBranchAddress("isoTrack_n_AllHits_TRACK", &isoTrack_n_AllHits_TRACK, &b_isoTrack_n_AllHits_TRACK);
   fChain->SetBranchAddress("isoTrack_n_AllHits_MISSING_INNER", &isoTrack_n_AllHits_MISSING_INNER, &b_isoTrack_n_AllHits_MISSING_INNER);
   fChain->SetBranchAddress("isoTrack_n_AllHits_MISSING_OUTER", &isoTrack_n_AllHits_MISSING_OUTER, &b_isoTrack_n_AllHits_MISSING_OUTER);
   fChain->SetBranchAddress("isoTrack_n_LostHits_TRACK", &isoTrack_n_LostHits_TRACK, &b_isoTrack_n_LostHits_TRACK);
   fChain->SetBranchAddress("isoTrack_n_LostHits_MISSING_INNER", &isoTrack_n_LostHits_MISSING_INNER, &b_isoTrack_n_LostHits_MISSING_INNER);
   fChain->SetBranchAddress("isoTrack_n_LostHits_MISSING_OUTER", &isoTrack_n_LostHits_MISSING_OUTER, &b_isoTrack_n_LostHits_MISSING_OUTER);
   fChain->SetBranchAddress("isoTrack_n_ValidPixelHits", &isoTrack_n_ValidPixelHits, &b_isoTrack_n_ValidPixelHits);
   fChain->SetBranchAddress("isoTrack_n_ValidStripHits", &isoTrack_n_ValidStripHits, &b_isoTrack_n_ValidStripHits);
   fChain->SetBranchAddress("isoTrack_n_LostPixelHits_TRACK", &isoTrack_n_LostPixelHits_TRACK, &b_isoTrack_n_LostPixelHits_TRACK);
   fChain->SetBranchAddress("isoTrack_n_LostPixelHits_MISSING_INNER", &isoTrack_n_LostPixelHits_MISSING_INNER, &b_isoTrack_n_LostPixelHits_MISSING_INNER);
   fChain->SetBranchAddress("isoTrack_n_LostPixelHits_MISSING_OUTER", &isoTrack_n_LostPixelHits_MISSING_OUTER, &b_isoTrack_n_LostPixelHits_MISSING_OUTER);
   fChain->SetBranchAddress("isoTrack_n_LostStripHits_TRACK", &isoTrack_n_LostStripHits_TRACK, &b_isoTrack_n_LostStripHits_TRACK);
   fChain->SetBranchAddress("isoTrack_n_LostStripHits_MISSING_INNER", &isoTrack_n_LostStripHits_MISSING_INNER, &b_isoTrack_n_LostStripHits_MISSING_INNER);
   fChain->SetBranchAddress("isoTrack_n_LostStripHits_MISSING_OUTER", &isoTrack_n_LostStripHits_MISSING_OUTER, &b_isoTrack_n_LostStripHits_MISSING_OUTER);
   Notify();
}

Bool_t MyTauClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyTauClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyTauClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyTauClass_cxx
