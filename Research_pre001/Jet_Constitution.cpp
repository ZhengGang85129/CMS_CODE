#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TF1.h"
#include "vector"
using namespace std;
const Int_t kMax = 20000;

const double hgcal_lower_eta = 1.5;
const double hgcal_upper_eta = 3.;

const Double_t R_cut = .4;
const Double_t constraint_coef_pt = .5;

const double canvas_width = 1200;
const double font_size = .035;
const double font_offset= .6;

const float Left_Margin = .12;
const float Right_Margin = .1;
const float Bottom_Margin = .12;
const float Top_Margin = .1;

const float marker_size = .8;
void Find_min(Double_t [] , Int_t , Int_t , Double_t & , int & ,bool &);
void Jet_Constitution()
{
//    gStyle->SetOptStat(0);
    gROOT->Reset();
//    gPad->SetFixedAspectRatio();
    TFile *f = new TFile("/Users/chenzhenggang/Desktop/local_analysis/Ntuple.root");
    TTree *mytree = (TTree*)f->Get("Tree");
    TCanvas *c1=new TCanvas("c1","c1",canvas_width,canvas_width);
    TCanvas *c2=new TCanvas("c2","c2",canvas_width,canvas_width);
    c1->SetMargin(Left_Margin ,Right_Margin ,Bottom_Margin ,Top_Margin);
    c2->SetMargin(Left_Margin ,Right_Margin ,Bottom_Margin ,Top_Margin);
    c1->SetCanvasSize(600,600);
    c2->SetCanvasSize(600,600);
    //Branch for genparticle
    Int_t genparticle_track;
    Int_t genparticle_charge[kMax];
    Int_t genparticle_pdgId[kMax];
    Int_t genparticle_status[kMax];
    Double_t genparticle_px[kMax];
    Double_t genparticle_py[kMax];
    Double_t genparticle_pz[kMax];
    Double_t genparticle_pt[kMax];
    Double_t genparticle_E[kMax];
    Double_t genparticle_mass[kMax];
    Double_t genparticle_mt[kMax];
    Double_t genparticle_phi[kMax];
    Double_t genparticle_theta[kMax];
    Double_t genparticle_eta[kMax];
    Double_t genparticle_rapidity[kMax];
    
    mytree->SetBranchAddress("genparticle_track",&genparticle_track);
    mytree->SetBranchAddress("genparticle_charge",&genparticle_charge[genparticle_track]);
    mytree->SetBranchAddress("genparticle_pdgId",&genparticle_pdgId[genparticle_track]);
    mytree->SetBranchAddress("genparticle_status",&genparticle_status[genparticle_track]);
    mytree->SetBranchAddress("genparticle_px",&genparticle_px[genparticle_track]);
    mytree->SetBranchAddress("genparticle_py",&genparticle_py[genparticle_track]);
    mytree->SetBranchAddress("genparticle_pz",&genparticle_pz[genparticle_track]);
    mytree->SetBranchAddress("genparticle_pt",&genparticle_pt[genparticle_track]);
    mytree->SetBranchAddress("genparticle_mass",&genparticle_mass[genparticle_track]);
    mytree->SetBranchAddress("genparticle_mt",&genparticle_mt[genparticle_track]);
    mytree->SetBranchAddress("genparticle_phi",&genparticle_phi[genparticle_track]);
    mytree->SetBranchAddress("genparticle_theta",&genparticle_theta[genparticle_track]);
    mytree->SetBranchAddress("genparticle_eta",&genparticle_eta[genparticle_track]);
    mytree->SetBranchAddress("genparticle_rapidity",&genparticle_rapidity[genparticle_track]);
    mytree->SetBranchAddress("genparticle_E",&genparticle_E[genparticle_track]);

    //Branch for AK4PFJets  
    
    Int_t ak4pfjets_track;
    Int_t ak4pfjets_charge[kMax];
    Int_t ak4pfjets_pdgId[kMax];
    Double_t ak4pfjets_px[kMax];
    Double_t ak4pfjets_py[kMax];
    Double_t ak4pfjets_pz[kMax];
    Double_t ak4pfjets_pt[kMax];
    Double_t ak4pfjets_mass[kMax];
    Double_t ak4pfjets_mt[kMax];
    Double_t ak4pfjets_phi[kMax];
    Double_t ak4pfjets_theta[kMax];
    Double_t ak4pfjets_eta[kMax];
    Double_t ak4pfjets_rapidity[kMax];
    Double_t ak4pfjets_E[kMax];
    Int_t ak4pfjets_nConstituents[kMax];
    Int_t ak4pfjets_neutralHadronMultiplicity[kMax];
    Int_t ak4pfjets_chargedHadronMultiplicity[kMax];
    Int_t ak4pfjets_chargedMultiplicity[kMax];
    Int_t ak4pfjets_neutralMultiplicity[kMax];
    std::vector<Int_t> *ak4pfjets_tagging;
    TBranch *ak4pfjets;
    mytree->SetBranchAddress("ak4pfjets_track", &ak4pfjets_track );
    mytree->SetBranchAddress("ak4pfjets_charge", &ak4pfjets_charge[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_pdgId", ak4pfjets_pdgId);
    mytree->SetBranchAddress("ak4pfjets_px", &ak4pfjets_px[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_py", &ak4pfjets_py[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_pz", &ak4pfjets_pz[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_pt", &ak4pfjets_pt[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_mass", &ak4pfjets_mass[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_mt", &ak4pfjets_mt[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_phi", &ak4pfjets_phi[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_theta", &ak4pfjets_theta[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_eta", &ak4pfjets_eta[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_rapidity", &ak4pfjets_rapidity[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_E", &ak4pfjets_E[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_nConstituents", &ak4pfjets_nConstituents[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_neutralHadronMultiplicity", &ak4pfjets_neutralHadronMultiplicity[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_chargedHadronMultiplicity", &ak4pfjets_chargedHadronMultiplicity[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_chargedMultiplicity", &ak4pfjets_chargedMultiplicity[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_neutralMultiplicity", &ak4pfjets_neutralMultiplicity[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_tagging", &ak4pfjets_tagging);

    //Branch for HGCal 
    
    Int_t hgcalclst_track;
    Int_t hgcalclst_size[kMax];
    Double_t hgcalclst_x[kMax];
    Double_t hgcalclst_y[kMax];
    Double_t hgcalclst_z[kMax];
    Double_t hgcalclst_E[kMax];
    Double_t hgcalclst_eta[kMax];
    Double_t hgcalclst_phi[kMax];
    Int_t hgcalclst_layer[kMax];
    mytree->SetBranchAddress("hgcalclst_track", &hgcalclst_track ) ;
    mytree->SetBranchAddress("hgcalclst_x",hgcalclst_x) ;
    mytree->SetBranchAddress("hgcalclst_y",hgcalclst_y) ;
    mytree->SetBranchAddress("hgcalclst_z",hgcalclst_z) ;
    mytree->SetBranchAddress("hgcalclst_E",hgcalclst_E) ;
    mytree->SetBranchAddress("hgcalclst_eta",hgcalclst_eta) ;
    mytree->SetBranchAddress("hgcalclst_phi",hgcalclst_phi) ;
    mytree->SetBranchAddress("hgcalclst_layer",hgcalclst_layer);

    //Branch for PFCandidate
    std::vector<Int_t> *PFConstituents_charge;
    std::vector<Int_t> *PFConstituents_pdgId;
    std::vector<Double_t> *PFConstituents_p;
    std::vector<Double_t> *PFConstituents_px;
    std::vector<Double_t> *PFConstituents_py;
    std::vector<Double_t> *PFConstituents_pz;
    std::vector<Double_t> *PFConstituents_pt;
    std::vector<Double_t> *PFConstituents_E;
    std::vector<Double_t> *PFConstituents_Et;
    std::vector<Double_t> *PFConstituents_mass;
    std::vector<Double_t> *PFConstituents_mt;
    std::vector<Double_t> *PFConstituents_phi;
    std::vector<Double_t> *PFConstituents_theta;
    std::vector<Double_t> *PFConstituents_eta;
    std::vector<Double_t> *PFConstituents_hcalEnergy;
    std::vector<Double_t> *PFConstituents_ecalEnergy;
    std::vector<Double_t> *PFConstituents_rapidity;

    TBranch *PFConstituents;
    //mytree->SetBranchAddress("PFConstituents_charge", &PFConstituents_charge,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_charge", &PFConstituents_charge,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_pdgId", &PFConstituents_pdgId,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_p", &PFConstituents_p,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_px", &PFConstituents_px,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_py", &PFConstituents_py,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_pz", &PFConstituents_pz,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_pt", &PFConstituents_pt,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_E", &PFConstituents_E,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_Et", &PFConstituents_Et,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_mass", &PFConstituents_mass,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_mt", &PFConstituents_mt,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_phi", &PFConstituents_phi,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_theta", &PFConstituents_theta,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_eta", &PFConstituents_eta,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_rapidity", &PFConstituents_rapidity,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_hcalEnergy", &PFConstituents_hcalEnergy,&PFConstituents ) ;
    mytree->SetBranchAddress("PFConstituents_ecalEnergy", &PFConstituents_ecalEnergy,&PFConstituents ) ;

    const int hgcallayer_size = 35;    
    const int MAX_NCONSTITUENTS=200;
//For Quark
    TH1D *Quark_nConstituents = new TH1D("Number of Constituents for q"," Quark nConstituents " , 40 , 0, MAX_NCONSTITUENTS  ) ;
//    TH1D *Quark_chargedHadronMultiplicity = new TH1D("Charged Hadron Multiplicity for Quark","Charged Hadron Multiplicity for Quark", 16 , 0, 80 ) ;
//    TH1D *Quark_neutralHadronMultiplicity = new TH1D("Neutral Hadron Multiplicity for Quark","Neutral Hadron Multiplicity for Quark", 16 , 0, 80 ) ;
    TH1D *Quark_neutralMultiplicity = new TH1D("Neutral Multiplicity for Quark","Neutral Multiplicity for Quark", 40 , 0, MAX_NCONSTITUENTS ) ;
    TH1D *Quark_chargedMultiplicity = new TH1D("Charged Multiplicity for Quark","Charged Multiplicity for Quark", 40 , 0, MAX_NCONSTITUENTS ) ;

//For Gluon    
    TH1D *Gluon_nConstituents = new TH1D("Number of Constituents for g"," Gluon nConstituents " , 40 , 0, MAX_NCONSTITUENTS ) ;
//    TH1D *Gluon_chargedHadronMultiplicity = new TH1D("Charged Hadron Multiplicity for Gluon","Charged Hadron Multiplicity for Gluon", 16 , 0, 80 ) ;
//    TH1D *Gluon_neutralHadronMultiplicity = new TH1D("Neutral Hadron Multiplicity for Gluon","Neutral Hadron Multiplicity for Gluon", 16 , 0, 80 ) ;
    TH1D *Gluon_neutralMultiplicity = new TH1D("Neutral Multiplicity for Gluon","Neutral Multiplicity for Gluon", 40 , 0, MAX_NCONSTITUENTS ) ;
    TH1D *Gluon_chargedMultiplicity = new TH1D("Charged Multiplicity for Gluon","Charged Multiplicity for Gluon", 40 , 0, MAX_NCONSTITUENTS ) ;

/*Setting of objects*/

//For Quark   
    Quark_nConstituents->GetXaxis()->SetTitle("nConstituents");
    Quark_nConstituents->GetYaxis()->SetTitle("Events");
    Quark_nConstituents->GetXaxis()->SetTitleSize(font_size);
    Quark_nConstituents->GetYaxis()->SetTitleSize(font_size); 
    Quark_nConstituents->GetXaxis()->SetTitleOffset(font_offset);
    Quark_nConstituents->GetYaxis()->SetTitleOffset(font_offset);
    Quark_nConstituents->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_nConstituents->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_nConstituents->SetFillColor(kRed);
    Quark_nConstituents->SetLineColor(kBlack);
/*
    Quark_neutralHadronMultiplicity->GetXaxis()->SetTitle("Neutral Hadron Multiplicity");
    Quark_neutralHadronMultiplicity->GetYaxis()->SetTitle("Events");
    Quark_neutralHadronMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Quark_neutralHadronMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Quark_neutralHadronMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Quark_neutralHadronMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Quark_neutralHadronMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_neutralHadronMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_neutralHadronMultiplicity->SetFillColor(kRed);  
    Quark_neutralHadronMultiplicity->SetLineColor(kBlack);

    Quark_chargedHadronMultiplicity->GetXaxis()->SetTitle("Charged Hadron Multiplicity");
    Quark_chargedHadronMultiplicity->GetYaxis()->SetTitle("Events");
    Quark_chargedHadronMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Quark_chargedHadronMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Quark_chargedHadronMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Quark_chargedHadronMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Quark_chargedHadronMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_chargedHadronMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_chargedHadronMultiplicity->SetFillColor(kRed);
    Quark_chargedHadronMultiplicity->SetLineColor(kBlack);
*/
    Quark_chargedMultiplicity->GetXaxis()->SetTitle("Charged Multiplicity");
    Quark_chargedMultiplicity->GetYaxis()->SetTitle("Events");
    Quark_chargedMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Quark_chargedMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Quark_chargedMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Quark_chargedMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Quark_chargedMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_chargedMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_chargedMultiplicity->SetFillColor(kRed);
    Quark_chargedMultiplicity->SetLineColor(kBlack);

    Quark_neutralMultiplicity->GetXaxis()->SetTitle("Neutral Multiplicity");
    Quark_neutralMultiplicity->GetYaxis()->SetTitle("Events");
    Quark_neutralMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Quark_neutralMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Quark_neutralMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Quark_neutralMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Quark_neutralMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_neutralMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_neutralMultiplicity->SetFillColor(kRed);
    Quark_neutralMultiplicity->SetLineColor(kBlack);


//For Gluon

    Gluon_nConstituents->GetXaxis()->SetTitle("nConstituents");
    Gluon_nConstituents->GetYaxis()->SetTitle("Events");
    Gluon_nConstituents->GetXaxis()->SetTitleSize(font_size);
    Gluon_nConstituents->GetYaxis()->SetTitleSize(font_size); 
    Gluon_nConstituents->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_nConstituents->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_nConstituents->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_nConstituents->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_nConstituents->SetFillColor(kGreen);
    Gluon_nConstituents->SetLineColor(kBlack);
/*
    Gluon_neutralHadronMultiplicity->GetXaxis()->SetTitle("Neutral Hadron Multiplicity");
    Gluon_neutralHadronMultiplicity->GetYaxis()->SetTitle("Events");
    Gluon_neutralHadronMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Gluon_neutralHadronMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Gluon_neutralHadronMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_neutralHadronMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_neutralHadronMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_neutralHadronMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_neutralHadronMultiplicity->SetFillColor(kGreen); 
    Gluon_neutralHadronMultiplicity->SetLineColor(kBlack);

    Gluon_chargedHadronMultiplicity->GetXaxis()->SetTitle("Charged Hadron Multiplicity");
    Gluon_chargedHadronMultiplicity->GetYaxis()->SetTitle("Events");
    Gluon_chargedHadronMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Gluon_chargedHadronMultiplicity->GetYaxis()->SetTitleSize(font_size);
    Gluon_chargedHadronMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_chargedHadronMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_chargedHadronMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_chargedHadronMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_chargedHadronMultiplicity->SetFillColor(kGreen); 
    Gluon_chargedHadronMultiplicity->SetLineColor(kBlack);
*/
    Gluon_chargedMultiplicity->GetXaxis()->SetTitle("Charged Multiplicity");
    Gluon_chargedMultiplicity->GetYaxis()->SetTitle("Events");
    Gluon_chargedMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Gluon_chargedMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Gluon_chargedMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_chargedMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_chargedMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_chargedMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_chargedMultiplicity->SetFillColor(kGreen);
    Gluon_chargedMultiplicity->SetLineColor(kBlack);

    Gluon_neutralMultiplicity->GetXaxis()->SetTitle("Neutral Multiplicity");
    Gluon_neutralMultiplicity->GetYaxis()->SetTitle("Events");
    Gluon_neutralMultiplicity->GetXaxis()->SetTitleSize(font_size);
    Gluon_neutralMultiplicity->GetYaxis()->SetTitleSize(font_size); 
    Gluon_neutralMultiplicity->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_neutralMultiplicity->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_neutralMultiplicity->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_neutralMultiplicity->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_neutralMultiplicity->SetFillColor(kGreen);
    Gluon_neutralMultiplicity->SetLineColor(kBlack);

    /*Setting of Histogram*/
    
    /*Setting for Selection*/
    Int_t nentries = (Int_t)mytree->GetEntries();
    
    for( Int_t ev = 0 ; ev < nentries ; ev++)
    {
      mytree->GetEntry(ev);
      for(int i = 0; i < genparticle_track  ; i ++)
      {
          if(TMath::Abs(genparticle_pdgId[ i ])< 4 && TMath::Abs(genparticle_pdgId[ i ]) > 0 && genparticle_status[ i ]  ==23  )
       {
           Double_t dR[ ak4pfjets_track ] ;
           Double_t minimum_DeltaR;
           int min_pos=0;
           bool flag=false;
           for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                  {
                      dR[ j ] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[ j ] - genparticle_eta[ i ]),2) + TMath::Power((ak4pfjets_phi[ j ] - genparticle_phi[ i ]), 2 ) );
                  }
           
           Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos , flag );
           if( minimum_DeltaR < 0.1 \
                   && ak4pfjets_pt[  min_pos ] > genparticle_pt[ min_pos ] * constraint_coef_pt\
                   && ak4pfjets_pt[  min_pos ] )
           {
               Double_t dR_; 
               for(Int_t k = 0 ; k < hgcalclst_track ; k ++ )
                {
                    dR_ = TMath::Sqrt(TMath::Power(hgcalclst_eta[ k ] - ak4pfjets_eta[ min_pos ] , 2 ) + TMath::Power( hgcalclst_phi[ k ] - ak4pfjets_phi[ min_pos ] , 2 ) );
                    if( dR_ < R_cut )
                    {
                    }
                }
               int charge_ctr=0; 
               int neutral_ctr=0;
               if(min_pos >0 )
                {
                    Int_t mindex = ak4pfjets_tagging->at(min_pos-1);
                    Int_t mindex_end = ak4pfjets_tagging->at(min_pos);
                    for(;mindex < mindex_end;mindex++)
                    {
                        if(PFConstituents_charge->at(mindex) !=0)
                        {
                            charge_ctr+=1;
                        }
                        else if(PFConstituents_charge->at(mindex) ==0 \
                                && PFConstituents_pt->at(mindex) >1 )
                        {
                            neutral_ctr+=1;
                        }
                    }
                }
               else
               {
                   for(int mindex=0;mindex<ak4pfjets_tagging->at(min_pos);mindex++)
                   {
                        if(PFConstituents_charge->at(mindex) !=0)
                        {
                            charge_ctr+=1;
                        }
                        else if(PFConstituents_charge->at(mindex) ==0 \
                                && PFConstituents_pt->at(mindex) >1 
                                )
                        {
                            neutral_ctr+=1;
                        }
                   }
               }
                Quark_nConstituents->Fill(neutral_ctr+charge_ctr);
                Quark_neutralMultiplicity->Fill(neutral_ctr);
                Quark_chargedMultiplicity->Fill(charge_ctr);
                }
       }
       else if(genparticle_pdgId[ i ] == 21 && genparticle_status[ i ] == 23 )
       {
           Double_t dR[ ak4pfjets_track ] ;
           Double_t minimum_DeltaR ;
           int min_pos=0;
           bool flag=false;

           for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                  {
                      dR[ j ] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[ j ] - genparticle_eta[ i ]),2) + TMath::Power((ak4pfjets_phi[ j ] - genparticle_phi[ i ]), 2 ) );
                  }
           Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos ,flag);
           if( minimum_DeltaR < 0.1   \
                   && ak4pfjets_pt[  min_pos ] > genparticle_pt[ min_pos ] * constraint_coef_pt\
                   && ak4pfjets_pt[  min_pos ] )
           {
               Double_t dR_ ;
               for(Int_t k = 0 ; k < hgcalclst_track ; k ++)
               {
                   dR_ = TMath::Sqrt(TMath::Power(hgcalclst_eta[ k ] - ak4pfjets_eta[ min_pos ] , 2 ) + TMath::Power( hgcalclst_phi[ k ] - ak4pfjets_phi[ min_pos ] , 2 ) );
                   if( dR_ < R_cut)
                   {
                   }
               }
               int charge_ctr=0; 
               int neutral_ctr=0;
               if(min_pos >0 )
                {
                    printf("%d\n",min_pos);
                    Int_t mindex = ak4pfjets_tagging->at(min_pos-1);
                    Int_t mindex_end = ak4pfjets_tagging->at(min_pos);
                    for(;mindex < mindex_end;mindex++)
                    {
                        if(PFConstituents_charge->at(mindex) !=0)
                        {
                            charge_ctr+=1;
                        }
                        else if(PFConstituents_charge->at(mindex) ==0 \
                                && PFConstituents_pt->at(mindex) >1 )
                        {
                            neutral_ctr+=1;
                        }
                    }
                }
               else
               {
                   for(int mindex=0;mindex<ak4pfjets_tagging->at(min_pos);mindex++)
                   {
                        if(PFConstituents_charge->at(mindex) !=0)
                        {
                            charge_ctr+=1;
                        }
                        else if(PFConstituents_charge->at(mindex) ==0 \
                                && PFConstituents_pt->at(mindex) >1 )
                        {
                            neutral_ctr+=1;
                        }
                   }
               }
                Gluon_nConstituents->Fill(neutral_ctr+charge_ctr);
                Gluon_neutralMultiplicity->Fill(neutral_ctr);
                Gluon_chargedMultiplicity->Fill(charge_ctr);
            }
       }
      }
    }
    c1->cd();
    Quark_nConstituents->Draw();
    c1->SaveAs("nConstituents_quark.pdf");
    c1->Clear();
/*
    Quark_chargedHadronMultiplicity->Draw();
    c1->SaveAs("Quark_chargedHadronMultiplicity.pdf");
    c1->Clear();
    Quark_neutralHadronMultiplicity->Draw();
    c1->SaveAs("Quark_neutralHadronMultiplicity.pdf");
    c1->Clear();
*/
    Quark_neutralMultiplicity->Draw();
    c1->SaveAs("Quark_neutralMultiplicity.pdf");
    c1->Clear();
    Quark_chargedMultiplicity->Draw();
    c1->SaveAs("Quark_chargedMultiplicity.pdf");
    c1->Clear();
    
    c2->cd();
    Gluon_nConstituents->Draw();
    c2->SaveAs("nConstituents_gluon.pdf");
    c2->Clear();
/*
    Gluon_chargedHadronMultiplicity->Draw();
    c2->SaveAs("Gluon_chargedHadronMultiplicity.pdf");
    c2->Clear();
    Gluon_neutralHadronMultiplicity->Draw();
    c2->SaveAs("Gluon_neutralHadronMultiplicity.pdf");
    c2->Clear();
*/
    Gluon_neutralMultiplicity->Draw();
    c2->SaveAs("Gluon_neutralMultiplicity.pdf");
    c2->Clear();
    Gluon_chargedMultiplicity->Draw();
    c2->SaveAs("Gluon_chargedMultiplicity.pdf");
    c2->Clear();
    
    c1->Close();
    c2->Close();
}
void Find_min(Double_t arr[] , const Int_t size , Int_t pos ,Double_t &min , int &i , bool &dsc)
{
    if(!dsc){
       min = arr[ 0 ];
       dsc = true;
    }
    if(pos != size - 1)
    {
        Find_min( arr, size , pos + 1 , min , i , dsc);
    }
    if( min > arr[ pos ])
    {
        min = arr[ pos ];
        i = pos ;
    }

}

