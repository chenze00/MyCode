#define MyTauClass_cxx
#include "MyCode/TupleProducer/interface/MyTauClass.h"
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <iostream>
#include <string.h>
//#include "plot1D.C"
#include "MyCode/TupleProducer/interface/GenLepton.h"
//#include <Math/LorentzVector.h>
#include <TColor.h>
#include <chrono>
#include <algorithm>

using namespace std;

float dPhi(float phi1,float phi2);// |phi1 - phi2| in 2pi range
bool sortpt( const vector<float_t>& v1,const vector<float_t>& v2 );//2d list sorted by row/column
void transpose(vector<vector<float_t> > &b);//matrix transpose

void MyTauClass::Loop(TString outputpath="")
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();//estimated time 375.327[s]
  // made with tau->MakeClass("MyTauClass")
  // to run it: 
  // .L MyTauClass.C++
  // MyTauClass aa
  // aa.Loop()
//   In a ROOT session, you can do:
//      root> .L MyTauClass.C
//      root> MyTauClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry{

  //notice that from now on there's an event preselection (i.e. genLepton_kind==5) which picks up 14% of the orginal data

   if (fChain == 0) return;
   Int_t ievent=0;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
// Book here my histograms
   
   //Ntuple output dataforML_(Tbranch)_(pfCand cone size)_(pt ratio)_(cone seed).root
    TString outputname=inputName+TString("__dataforML_flattend2_dR02_jetpt_jet.root");
    if(outputpath.Length()>0 && !outputpath.EndsWith("/")) outputpath.Append("/");
    TFile ML_file(outputpath+outputname,"RECREATE");
    TTree ML_data("ML_data", "Example N-Tuple");
   
    int dm_gen, dm_HPS;
    float_t tempSumW;
   //vector<float> ptratio_pfCand;
   // ML_data.Branch("pt_pfCand_ratio",&ptratio_pfCand);
    
    int n_photons,n_chargedHadrons,n_neutralHadrons,n_electrons,n_muons;
    float pt_ratio_Lphoton,pt_ratio_NLphoton,pt_ratio_NNLphoton,pt_ratio_NNNLphoton,deltaEta_Lphoton,deltaEta_NLphoton,deltaEta_NNLphoton,deltaEta_NNNLphoton,deltaPhi_Lphoton,deltaPhi_NLphoton,deltaPhi_NNLphoton,deltaPhi_NNNLphoton;
    float pt_ratio_LchargedHadron,pt_ratio_NLchargedHadron,pt_ratio_NNLchargedHadron,deltaEta_LchargedHadron,deltaEta_NLchargedHadron,deltaEta_NNLchargedHadron,deltaPhi_LchargedHadron,deltaPhi_NLchargedHadron,deltaPhi_NNLchargedHadron;
    float pt_ratio_LneutralHadron,pt_ratio_NLneutralHadron,deltaEta_LneutralHadron,deltaEta_NLneutralHadron,deltaPhi_LneutralHadron,deltaPhi_NLneutralHadron;
    float pt_ratio_Lelectron,pt_ratio_NLelectron,pt_ratio_NNLelectron,pt_ratio_NNNLelectron,deltaEta_Lelectron,deltaEta_NLelectron,deltaEta_NNLelectron,deltaEta_NNNLelectron,deltaPhi_Lelectron,deltaPhi_NLelectron,deltaPhi_NNLelectron,deltaPhi_NNNLelectron;
    float pt_ratio_Lmuon,pt_ratio_NLmuon,deltaEta_Lmuon,deltaEta_NLmuon,deltaPhi_Lmuon,deltaPhi_NLmuon;
    
    ML_data.Branch("dm_gen",&dm_gen,"dm_gen/I");
    ML_data.Branch("dm_HPS",&dm_HPS,"dm_HPS/I");
    ML_data.Branch("dR_weighted_by_ptratio",&tempSumW,"dR_weighted_by_ptratio/F");

    ML_data.Branch("n_photons",&n_photons,"n_photons/I");
    ML_data.Branch("pt_ratio_Lphoton",&pt_ratio_Lphoton,"pt_ratio_Lphoton/F");
    ML_data.Branch("pt_ratio_NLphoton",&pt_ratio_NLphoton,"pt_ratio_NLphoton/F");
    ML_data.Branch("pt_ratio_NNLphoton",&pt_ratio_NNLphoton,"pt_ratio_NNLphoton/F");
    ML_data.Branch("pt_ratio_NNNLphoton",&pt_ratio_NNNLphoton,"pt_ratio_NNNLphoton/F");
    ML_data.Branch("deltaEta_Lphoton",&deltaEta_Lphoton,"deltaEta_Lphoton/F");
    ML_data.Branch("deltaEta_NLphoton",&deltaEta_NLphoton,"deltaEta_NLphoton/F");
    ML_data.Branch("deltaEta_NNLphoton",&deltaEta_NNLphoton,"deltaEta_NNLphoton/F");
    ML_data.Branch("deltaEta_NNNLphoton",&deltaEta_NNNLphoton,"deltaEta_NNNLphoton/F");    
    ML_data.Branch("deltaPhi_Lphoton",&deltaPhi_Lphoton,"deltaPhi_Lphoton/F");
    ML_data.Branch("deltaPhi_NLphoton",&deltaPhi_NLphoton,"deltaPhi_NLphoton/F");
    ML_data.Branch("deltaPhi_NNLphoton",&deltaPhi_NNLphoton,"deltaPhi_NNLphoton/F");
    ML_data.Branch("deltaPhi_NNNLphoton",&deltaPhi_NNNLphoton,"deltaPhi_NNNLphoton/F");
    //ML_data.Branch("deltaR_Lphoton",&deltaR_Lphoton,"deltaR_Lphoton/F");
    //ML_data.Branch("deltaR_NLphoton",&deltaR_NLphoton,"deltaR_NLphoton/F");

    ML_data.Branch("n_chargedHadrons",&n_chargedHadrons,"n_chargedHadrons/I");
    ML_data.Branch("pt_ratio_LchargedHadron",&pt_ratio_LchargedHadron,"pt_ratio_LchargedHadron/F");
    ML_data.Branch("pt_ratio_NLchargedHadron",&pt_ratio_NLchargedHadron,"pt_ratio_NLchargedHadron/F");
    ML_data.Branch("pt_ratio_NNLchargedHadron",&pt_ratio_NNLchargedHadron,"pt_ratio_NNLchargedHadron/F");
    ML_data.Branch("deltaEta_LchargedHadron",&deltaEta_LchargedHadron,"deltaEta_LchargedHadron/F");
    ML_data.Branch("deltaEta_NLchargedHadron",&deltaEta_NLchargedHadron,"deltaEta_NLchargedHadron/F");
    ML_data.Branch("deltaEta_NNLchargedHadron",&deltaEta_NNLchargedHadron,"deltaEta_NNLchargedHadron/F");
    ML_data.Branch("deltaPhi_LchargedHadron",&deltaPhi_LchargedHadron,"deltaPhi_LchargedHadron/F");
    ML_data.Branch("deltaPhi_NLchargedHadron",&deltaPhi_NLchargedHadron,"deltaPhi_NLchargedHadron/F");
    ML_data.Branch("deltaPhi_NNLchargedHadron",&deltaPhi_NNLchargedHadron,"deltaPhi_NNLchargedHadron/F");
    //ML_data.Branch("deltaR_LchargedHadron",&deltaR_LchargedHadron,"deltaR_LchargedHadron/F");    
    //ML_data.Branch("deltaR_NLchargedHadron",&deltaR_NLchargedHadron,"deltaR_NLchargedHadron/F");    

    ML_data.Branch("n_neutralHadrons",&n_neutralHadrons,"n_neutralHadrons/I");
    ML_data.Branch("pt_ratio_LneutralHadron",&pt_ratio_LneutralHadron,"pt_ratio_LneutralHadron/F");
    ML_data.Branch("pt_ratio_NLneutralHadron",&pt_ratio_NLneutralHadron,"pt_ratio_NLneutralHadron/F");
    ML_data.Branch("deltaEta_LneutralHadron",&deltaEta_LneutralHadron,"deltaEta_LneutralHadron/F");
    ML_data.Branch("deltaEta_NLneutralHadron",&deltaEta_NLneutralHadron,"deltaEta_NLneutralHadron/F");
    ML_data.Branch("deltaPhi_LneutralHadron",&deltaPhi_LneutralHadron,"deltaPhi_LneutralHadron/F");
    ML_data.Branch("deltaPhi_NLneutralHadron",&deltaPhi_NLneutralHadron,"deltaPhi_NLneutralHadron/F");
    //ML_data.Branch("deltaR_LneutralHadron",&deltaR_LneutralHadron,"deltaR_LneutralHadron/F");    
    //ML_data.Branch("deltaR_NLneutralHadron",&deltaR_NLneutralHadron,"deltaR_NLneutralHadron/F");

    ML_data.Branch("n_electrons",&n_electrons,"n_electrons/I");
    ML_data.Branch("pt_ratio_Lelectron",&pt_ratio_Lelectron,"pt_ratio_Lelectron/F");
    ML_data.Branch("pt_ratio_NLelectron",&pt_ratio_NLelectron,"pt_ratio_NLelectron/F");
    ML_data.Branch("pt_ratio_NNLelectron",&pt_ratio_NNLelectron,"pt_ratio_NNLelectron/F");
    ML_data.Branch("pt_ratio_NNNLelectron",&pt_ratio_NNNLelectron,"pt_ratio_NNNLelectron/F");
    ML_data.Branch("deltaEta_Lelectron",&deltaEta_Lelectron,"deltaEta_Lelectron/F");
    ML_data.Branch("deltaEta_NLelectron",&deltaEta_NLelectron,"deltaEta_NLelectron/F");
    ML_data.Branch("deltaEta_NNLelectron",&deltaEta_NNLelectron,"deltaEta_NNLelectron/F");
    ML_data.Branch("deltaEta_NNNLelectron",&deltaEta_NNNLelectron,"deltaEta_NNNLelectron/F");
    ML_data.Branch("deltaPhi_Lelectron",&deltaPhi_Lelectron,"deltaPhi_Lelectron/F");
    ML_data.Branch("deltaPhi_NLelectron",&deltaPhi_NLelectron,"deltaPhi_NLelectron/F");  
    ML_data.Branch("deltaPhi_NNLelectron",&deltaPhi_NNLelectron,"deltaPhi_NNLelectron/F");
    ML_data.Branch("deltaPhi_NNNLelectron",&deltaPhi_NNNLelectron,"deltaPhi_NNNLelectron/F"); 
    //ML_data.Branch("deltaR_Lelectron",&deltaR_Lelectron,"deltaR_Lelectron/F");
    //ML_data.Branch("deltaR_NLelectron",&deltaR_NLelectron,"deltaR_NLelectron/F");

    ML_data.Branch("n_muons",&n_muons,"n_muons/I");
    ML_data.Branch("pt_ratio_Lmuon",&pt_ratio_Lmuon,"pt_ratio_Lmuon/F");
    ML_data.Branch("pt_ratio_NLmuon",&pt_ratio_NLmuon,"pt_ratio_NLmuon/F");
    ML_data.Branch("deltaEta_Lmuon",&deltaEta_Lmuon,"deltaEta_Lmuon/F");
    ML_data.Branch("deltaEta_NLmuon",&deltaEta_NLmuon,"deltaEta_NLmuon/F");
    ML_data.Branch("deltaPhi_Lmuon",&deltaPhi_Lmuon,"deltaPhi_Lmuon/F");
    ML_data.Branch("deltaPhi_NLmuon",&deltaPhi_NLmuon,"deltaPhi_NLmuon/F");
    //ML_data.Branch("deltaR_Lmuon",&deltaR_Lmuon,"deltaR_Lmuon/F");
    //ML_data.Branch("deltaR_NLmuon",&deltaR_NLmuon,"deltaR_NLmuon/F");


    TH2F *dR_vs_gentau_pt=new TH2F("dR_vs_gentau_pt","dR(gen_tau,reco_jet) vs gentau_pt",50,0.,250.,50,0.,.5);//
    TH2F *dR_vs_gentau_pt_dmgen0=new TH2F("dR_vs_gentau_pt_dmgen0","dR(gen_tau,reco_jet) vs gentau_pt in gen DM 0",50,0.,250.,50,0.,.5);
    TH2F *dR_vs_gentau_pt_dmgen1=new TH2F("dR_vs_gentau_pt_dmgen1","dR(gen_tau,reco_jet) vs gentau_pt in gen DM 1",50,0.,250.,50,0.,.5);
    TH2F *dR_vs_gentau_pt_dmgen2=new TH2F("dR_vs_gentau_pt_dmgen2","dR(gen_tau,reco_jet) vs gentau_pt in gen DM 2",50,0.,250.,50,0.,.5);
    TH2F *dR_vs_gentau_pt_dmgen10=new TH2F("dR_vs_gentau_pt_dmgen10","dR(gen_tau,reco_jet) vs gentau_pt in gen DM 10",50,0.,250.,50,0.,.5);
    TH2F *dR_vs_gentau_pt_dmgen11=new TH2F("dR_vs_gentau_pt_dmgen11","dR(gen_tau,reco_jet) vs gentau_pt in gen DM 11",50,0.,250.,50,0.,.5);

    TH2F *dR1_vs_gentau_pt=new TH2F("dR1_vs_gentau_pt","dR(gen_tau,LchargedHadron) vs gentau_pt",50,0.,250.,50,0.,.5);//
    TH2F *dR1_vs_gentau_pt_dmgen0=new TH2F("dR1_vs_gentau_pt_dmgen0","dR(gen_tau,LchargedHadron) vs gentau_pt in gen DM 0",50,0.,250.,50,0.,.5);
    TH2F *dR1_vs_gentau_pt_dmgen1=new TH2F("dR1_vs_gentau_pt_dmgen1","dR(gen_tau,LchargedHadron) vs gentau_pt in gen DM 1",50,0.,250.,50,0.,.5);
    TH2F *dR1_vs_gentau_pt_dmgen2=new TH2F("dR1_vs_gentau_pt_dmgen2","dR(gen_tau,LchargedHadron) vs gentau_pt in gen DM 2",50,0.,250.,50,0.,.5);
    TH2F *dR1_vs_gentau_pt_dmgen10=new TH2F("dR1_vs_gentau_pt_dmgen10","dR(gen_tau,LchargedHadron) vs gentau_pt in gen DM 10",50,0.,250.,50,0.,.5);
    TH2F *dR1_vs_gentau_pt_dmgen11=new TH2F("dR1_vs_gentau_pt_dmgen11","dR(gen_tau,LchargedHadron) vs gentau_pt in gen DM 11",50,0.,250.,50,0.,.5);
    
    TH1F *visibleP4_pt=new TH1F("visibleP4_pt","visiblePt of genLepton",100,0.,100.);
    TH1F *radiatedP4_pt=new TH1F("radiatedP4_pt","radiatedPt of genLepton",100,0.,100.);
    
    TH1F *pt_genLep = new TH1F("pt_genLep_and_pt_tau","visible pt of  generated tau",100,0.,100.);
    TH1F *pt_HPSTau = new TH1F("pt_HPSTau","pt of  HPS reconstructed tau",100,0.,100.);
    TH1F *pt_jet = new TH1F("pt_jet","pt of jet",100,0.,100.);
  
// start of the loop
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      ievent++;
      if(ievent%30000==0) cout << "processed " << float(ievent)/nentries*100.  <<"%"<< endl;
      //cout<<"ievent="<<ievent<<std::endl;
      int dm=tau_decayMode; // get decay mode
       dm_HPS=dm;
      Int_t leng = pfCand_pt->size();//number of pfCand

      if(leng<=0)continue;
      if(genLepton_kind!=5)continue;
      if(jet_pt<=0)continue;
      if(tau_pt<=0)continue;

	//event selection    && genLepton_kind==5 && tau_pt>0
	
	auto genLeptons = reco_tau::gen_truth::GenLepton::fromRootTuple(//<std::vector>(
										     genLepton_lastMotherIndex,
										     *genParticle_pdgId,
										     *genParticle_mother,
										     *genParticle_charge,
										     *genParticle_isFirstCopy,
										     *genParticle_isLastCopy,
										     *genParticle_pt,
										     *genParticle_eta,
										     *genParticle_phi,
										     *genParticle_mass,
										     *genParticle_vtx_x,
										     *genParticle_vtx_y,
										     *genParticle_vtx_z);
	
	int nChargedHadrons = genLeptons.nChargedHadrons();
	int nNeutralHadrons = genLeptons.nNeutralHadrons();       
	dm_gen = (nChargedHadrons-1)*5 + nNeutralHadrons;
       
       pt_genLep->Fill(genLepton_vis_pt);
       pt_HPSTau->Fill(tau_pt);
       pt_jet->Fill(jet_pt);
	
	Float_t deltaR = 0, deltaphi = 0, deltaeta = 0;//deltaR between pfCand and genTau	  
	vector<float>pt_photons,pt_chargedHadrons,pt_neutralHadrons,pt_electrons,pt_muons;//actually it's ptratio :)
	vector<float>deltaEta_photons,deltaEta_chargedHadrons,deltaEta_neutralHadrons,deltaEta_electrons,deltaEta_muons;
	vector<float>deltaPhi_photons,deltaPhi_chargedHadrons,deltaPhi_neutralHadrons,deltaPhi_electrons,deltaPhi_muons;
	n_photons=0,n_chargedHadrons=0,n_neutralHadrons=0,n_electrons=0,n_muons=0;
	float dRSumWeighted=0,dRSumWeightedId1=0,dRSumWeightedId2=0,dRSumWeightedId3=0,dRSumWeightedId4=0,dRSumWeightedId5=0;

	Float_t Deta_genTau_jet=abs(genLepton_vis_eta - jet_eta);
	Float_t Dphi_genTau_jet=abs(dPhi(genLepton_vis_phi, jet_phi));
	Float_t dR_genTau_jet=sqrt(pow(Deta_genTau_jet,2) + pow(Dphi_genTau_jet,2));
	dR_vs_gentau_pt->Fill(genLepton_vis_pt,dR_genTau_jet);
   if(dm_gen==0)dR_vs_gentau_pt_dmgen0->Fill(genLepton_vis_pt,dR_genTau_jet);
   if(dm_gen==1)dR_vs_gentau_pt_dmgen1->Fill(genLepton_vis_pt,dR_genTau_jet );
   if(dm_gen==2||dm_gen==3||dm_gen==4||dm_gen==5)dR_vs_gentau_pt_dmgen2->Fill(genLepton_vis_pt,dR_genTau_jet );
   if(dm_gen==10)dR_vs_gentau_pt_dmgen10->Fill( genLepton_vis_pt,dR_genTau_jet );
   if(dm_gen==11)dR_vs_gentau_pt_dmgen11->Fill( genLepton_vis_pt,dR_genTau_jet );
	
	Float_t Lh_eta=-10,Lh_phi=-10,Lh_pt=-10;
	for(Int_t j=0; j<leng; j++){
	  if((*pfCand_particleType)[j]==1){
	    if((*pfCand_pt)[j]>Lh_pt){Lh_pt=(*pfCand_pt)[j];Lh_eta=(*pfCand_eta)[j];Lh_phi=(*pfCand_phi)[j];}
	  }
	}
	if(Lh_pt==-10){Lh_eta=jet_eta;Lh_phi=jet_phi;}

   if(Lh_pt!=-10){    
   Float_t Deta_genTau_Lh=abs(genLepton_vis_eta - Lh_eta);
	Float_t Dphi_genTau_Lh=abs(dPhi(genLepton_vis_phi, Lh_phi));
	Float_t dR_genTau_Lh=sqrt(pow(Deta_genTau_Lh,2) + pow(Dphi_genTau_Lh,2));
	dR1_vs_gentau_pt->Fill( genLepton_vis_pt,dR_genTau_Lh);
   if(dm_gen==0)dR1_vs_gentau_pt_dmgen0->Fill( genLepton_vis_pt,dR_genTau_Lh);
   if(dm_gen==1)dR1_vs_gentau_pt_dmgen1->Fill( genLepton_vis_pt,dR_genTau_Lh );
   if(dm_gen==2||dm_gen==3||dm_gen==4||dm_gen==5)dR1_vs_gentau_pt_dmgen2->Fill( genLepton_vis_pt,dR_genTau_Lh );
   if(dm_gen==10)dR1_vs_gentau_pt_dmgen10->Fill( genLepton_vis_pt,dR_genTau_Lh );
   if(dm_gen==11)dR1_vs_gentau_pt_dmgen11->Fill( genLepton_vis_pt,dR_genTau_Lh );
   }

           Float_t sumPT_pfCand=0;
        for(Int_t j=0; j<leng; j++){//run over pfCand in an event
	  int particleId=(*pfCand_particleType)[j];
            
     Float_t ConeSeed_eta=jet_eta, ConeSeed_phi=jet_phi;
	  deltaeta = abs(ConeSeed_eta - (*pfCand_eta)[j]);//genLepton_vis_eta OR jet_eta OR Lh_eta
	  deltaphi = abs(dPhi(ConeSeed_phi,(*pfCand_phi)[j]));//genLepton_vis_phi OR jet_phi OR Lh_eta
	  deltaR=sqrt(deltaeta*deltaeta + deltaphi*deltaphi);	    
        
	  if((*pfCand_pt)[j]>0. && deltaR<0.2){  //cone size 
	    if(particleId==1){pt_chargedHadrons.push_back((*pfCand_pt)[j]/jet_pt);deltaEta_chargedHadrons.push_back(deltaeta);deltaPhi_chargedHadrons.push_back(deltaphi);
	      dRSumWeightedId1=dRSumWeightedId1+(*pfCand_pt)[j]/jet_pt*deltaR;n_chargedHadrons++;}
	    if(particleId==2){pt_electrons.push_back((*pfCand_pt)[j]/jet_pt);deltaEta_electrons.push_back(deltaeta);deltaPhi_electrons.push_back(deltaphi);
	      dRSumWeightedId2=dRSumWeightedId2+(*pfCand_pt)[j]/jet_pt*deltaR;n_electrons++;}
	    if(particleId==3){pt_muons.push_back((*pfCand_pt)[j]/jet_pt);deltaEta_muons.push_back(deltaeta);deltaPhi_muons.push_back(deltaphi);
	      dRSumWeightedId3=dRSumWeightedId3+(*pfCand_pt)[j]/jet_pt*deltaR;n_muons++;}
	    if(particleId==4){pt_photons.push_back((*pfCand_pt)[j]/jet_pt);deltaEta_photons.push_back(deltaeta);deltaPhi_photons.push_back(deltaphi);
	      dRSumWeightedId4=dRSumWeightedId4+(*pfCand_pt)[j]/jet_pt*deltaR;n_photons++;}
	    if(particleId==5){pt_neutralHadrons.push_back((*pfCand_pt)[j]/jet_pt);deltaEta_neutralHadrons.push_back(deltaeta);deltaPhi_neutralHadrons.push_back(deltaphi);
	      dRSumWeightedId5=dRSumWeightedId5+(*pfCand_pt)[j]/jet_pt*deltaR;n_neutralHadrons++;}
	    dRSumWeighted=dRSumWeighted+(*pfCand_pt)[j]/jet_pt*deltaR;
	  }
	}//end pfCand loop

       vector<vector<float_t>> vect_chargedHadron={pt_chargedHadrons, deltaEta_chargedHadrons, deltaPhi_chargedHadrons};
       vector<vector<float_t>> vect_electron={pt_electrons, deltaEta_electrons, deltaPhi_electrons};
       vector<vector<float_t>> vect_muon={pt_muons, deltaEta_muons, deltaPhi_muons};
       vector<vector<float_t>> vect_photon={pt_photons, deltaEta_photons, deltaPhi_photons};
       vector<vector<float_t>> vect_neutralHadron={pt_neutralHadrons, deltaEta_neutralHadrons, deltaPhi_neutralHadrons};
       vector<float_t> vect_NaN={-10,-10,-10};
       transpose(vect_chargedHadron);for(int i=0;i<4;i++){vect_chargedHadron.push_back(vect_NaN);}
       transpose(vect_electron);for(int i=0;i<4;i++){vect_electron.push_back(vect_NaN);}
       transpose(vect_muon);for(int i=0;i<4;i++){vect_muon.push_back(vect_NaN);}
       transpose(vect_photon);for(int i=0;i<4;i++){vect_photon.push_back(vect_NaN);}
       transpose(vect_neutralHadron);for(int i=0;i<4;i++){vect_neutralHadron.push_back(vect_NaN);}
       
       sort(vect_chargedHadron.begin(),vect_chargedHadron.end(),sortpt);
       sort(vect_electron.begin(),vect_electron.end(),sortpt);
       sort(vect_muon.begin(),vect_muon.end(),sortpt);
       sort(vect_photon.begin(),vect_photon.end(),sortpt);
       sort(vect_neutralHadron.begin(),vect_neutralHadron.end(),sortpt);
       
       pt_ratio_LchargedHadron=vect_chargedHadron[0][0];
       pt_ratio_NLchargedHadron=vect_chargedHadron[1][0];
       pt_ratio_NNLchargedHadron=vect_chargedHadron[2][0];
       deltaEta_LchargedHadron=vect_chargedHadron[0][1];
       deltaEta_NLchargedHadron=vect_chargedHadron[1][1];
       deltaEta_NNLchargedHadron=vect_chargedHadron[2][1];
       deltaPhi_LchargedHadron=vect_chargedHadron[0][2];
       deltaPhi_NLchargedHadron=vect_chargedHadron[1][2];
       deltaPhi_NNLchargedHadron=vect_chargedHadron[2][2];

       pt_ratio_Lelectron=vect_electron[0][0];
       pt_ratio_NLelectron=vect_electron[1][0];
       pt_ratio_NNLelectron=vect_electron[2][0];
       pt_ratio_NNNLelectron=vect_electron[3][0];
       deltaEta_Lelectron=vect_electron[0][1];
       deltaEta_NLelectron=vect_electron[1][1];
       deltaEta_NNLelectron=vect_electron[2][1];
       deltaEta_NNNLelectron=vect_electron[3][1];
       deltaPhi_Lelectron=vect_electron[0][2];
       deltaPhi_NLelectron=vect_electron[1][2];
       deltaPhi_NNLelectron=vect_electron[2][2];
       deltaPhi_NNNLelectron=vect_electron[3][2];
       
       pt_ratio_Lmuon=vect_muon[0][0];
       pt_ratio_NLmuon=vect_muon[1][0];
       deltaEta_Lmuon=vect_muon[0][1];
       deltaEta_NLmuon=vect_muon[1][1];
       deltaPhi_Lmuon=vect_muon[0][2];
       deltaPhi_NLmuon=vect_muon[1][2];    

       pt_ratio_Lphoton=vect_photon[0][0];
       pt_ratio_NLphoton=vect_photon[1][0];
       pt_ratio_NNLphoton=vect_photon[2][0];
       pt_ratio_NNNLphoton=vect_photon[3][0];
       deltaEta_Lphoton=vect_photon[0][1];
       deltaEta_NLphoton=vect_photon[1][1];
       deltaEta_NNLphoton=vect_photon[2][1];
       deltaEta_NNNLphoton=vect_photon[3][1];
       deltaPhi_Lphoton=vect_photon[0][2];
       deltaPhi_NLphoton=vect_photon[1][2];
       deltaPhi_NNLphoton=vect_photon[2][2];
       deltaPhi_NNNLphoton=vect_photon[3][2];

       pt_ratio_LneutralHadron=vect_neutralHadron[0][0];
       pt_ratio_NLneutralHadron=vect_neutralHadron[1][0];
       deltaEta_LneutralHadron=vect_neutralHadron[0][1];
       deltaEta_NLneutralHadron=vect_neutralHadron[1][1];
       deltaPhi_LneutralHadron=vect_neutralHadron[0][2];
       deltaPhi_NLneutralHadron=vect_neutralHadron[1][2];
       //float_t tempSumWId1,tempSumWId2,tempSumWId3,tempSumWId4,tempSumWId5;
	tempSumW=dRSumWeighted;
	//tempSumWId1=dRSumWeightedId1;
	//tempSumWId2=dRSumWeightedId2;
	//tempSumWId3=dRSumWeightedId3;
	//tempSumWId4=dRSumWeightedId4;
	//tempSumWId5=dRSumWeightedId5;	
	
	ML_data.Fill();
      
   }   // end loop on events
   //write in file

   ML_data.Write();
   dR_vs_gentau_pt->Write();
   dR_vs_gentau_pt_dmgen0->Write();
   dR_vs_gentau_pt_dmgen1->Write();
   dR_vs_gentau_pt_dmgen2->Write();
   dR_vs_gentau_pt_dmgen10->Write();
   dR_vs_gentau_pt_dmgen11->Write();
   dR1_vs_gentau_pt->Write();
   dR1_vs_gentau_pt_dmgen0->Write();
   dR1_vs_gentau_pt_dmgen1->Write();
   dR1_vs_gentau_pt_dmgen2->Write();
   dR1_vs_gentau_pt_dmgen10->Write();
   dR1_vs_gentau_pt_dmgen11->Write();
    pt_genLep->Write();
    pt_HPSTau->Write();
    pt_jet->Write();
   ML_file.Close();

   //   TFile out_file("dataforplot.root","RECREATE");
   //  dR_vs_jet_pt->Write();
   // out_file.Close();
   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

   std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1e6 << "[s]" << std::endl;
    /*
    TCanvas *c1=new TCanvas("c1");
    THStack *stackplot = new THStack("stackplot","pt of visible generated tau and HPS reconstructed tau");
    TLegend *mylegend1 = new TLegend(0.6,0.6,0.8,0.8,"");
    pt_genLep->SetLineColor(kBlue);
    pt_HPSTau->SetLineColor(kOrange);
    pt_jet->SetLineColor(kRed);
    stackplot->Add(pt_genLep);
    stackplot->Add(pt_HPSTau);
    stackplot->Add(pt_jet);
    mylegend1->AddEntry(pt_genLep,"generated tau");
    mylegend1->AddEntry(pt_HPSTau,"HPS reconstructed tau");
    mylegend1->AddEntry(pt_jet,"jet");
    stackplot->Draw("nostack");
    mylegend1->Draw();
    c1->SaveAs("pt_gentau__and_HPS_tau_and_jet.pdf");
    */
    
}
float dPhi(float phi1, float phi2) {
    float result = phi1 - phi2;
    while (result > TMath::Pi()) result -= float(2.*TMath::Pi());
    while (result <= -TMath::Pi()) result += float(2.*TMath::Pi());
    float absresult=TMath::Abs(result);
    return absresult;
}

 bool sortpt( const vector<float_t>& v1,const vector<float_t>& v2 ) {
    return v1[0] > v2[0];
}

void transpose(vector<vector<float_t> > &b)
{
    if (b.size() == 0)
        return;

    vector<vector<float_t> > trans_vec(b[0].size(), vector<float_t>());

    for (unsigned i = 0; i < b.size(); i++)
    {
        for (unsigned j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }
    b = trans_vec;    // <--- reassign here
}
