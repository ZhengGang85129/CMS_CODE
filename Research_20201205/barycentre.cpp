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

using namespace std;
const Int_t kMax = 20000;
const int selected_q = 9;
const int selected_g = 2;

const double hgcal_lower_eta = 1.5;
const double hgcal_upper_eta = 3.0;

const double part_dR_cut =0.1;
const Double_t R_cut = 0.4;
const Double_t constraint_coef_pt = .5;

const int Number_of_HIST= 50;

const double canvas_width = 900;
const double font_size = 0.055;
const double font_offset= 0.9;

const float Left_Margin = 0.13;
const float Right_Margin = 0.1;
const float Bottom_Margin = 0.13;
const float Top_Margin = 0.1;

const float marker_size = 0.8;
const int hgcallayer_max=50;
void Find_min(Double_t [] , Int_t , Int_t , Double_t & , int & ,bool &);
void Normalization(Double_t [] , const Int_t);
const int Gluon_ID = 21;
const int Quark_lower_ID=1;
const int Quark_upper_ID=3;
const int Electron_ID=11;
const int Photon_ID=22;
double min_z = 322.102722;
double max_z = min_z+200;

void barycentre()
{
    delete gROOT->FindObject("Radiance_31");
    TSeqCollection* canvases = gROOT->GetListOfCanvases();
    TIter next(gROOT->GetListOfCanvases());
    while(c = (TCanvas*)next())
    {
        delete c;
    }
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
    mytree->SetBranchAddress("genparticle_pt",&genparticle_pt[genparticle_track]);
    mytree->SetBranchAddress("genparticle_phi",&genparticle_phi[genparticle_track]);
    mytree->SetBranchAddress("genparticle_eta",&genparticle_eta[genparticle_track]);
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
    std::vector<Int_t> *ak4pfjets_tagging;
    TBranch *ak4pfjets;
    mytree->SetBranchAddress("ak4pfjets_track", &ak4pfjets_track );
    mytree->SetBranchAddress("ak4pfjets_charge", &ak4pfjets_charge[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_pdgId", &ak4pfjets_pdgId );
    mytree->SetBranchAddress("ak4pfjets_pt", &ak4pfjets_pt[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_phi", &ak4pfjets_phi[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_eta", &ak4pfjets_eta[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_E", &ak4pfjets_E[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_tagging", &ak4pfjets_tagging);
    //Branch for HGCal 
    
    Int_t LayerCluster_track;
    Int_t LayerCluster_layer[kMax];
    Int_t LayerCluster_size[kMax];
    Double_t LayerCluster_x[kMax];
    Double_t LayerCluster_y[kMax];
    Double_t LayerCluster_z[kMax];
    Double_t LayerCluster_E[kMax];
    Double_t LayerCluster_eta[kMax];
    Double_t LayerCluster_phi[kMax];
    Double_t LayerCluster_mask_EM[kMax];
    Double_t LayerCluster_mask_HAD[kMax];
    Double_t LayerCluster_mask_Trk[kMax];
    Double_t LayerCluster_mask_MIP[kMax];

    mytree->SetBranchAddress("Cluster_track", &LayerCluster_track ) ;
    mytree->SetBranchAddress("Cluster_z",LayerCluster_z) ;
    mytree->SetBranchAddress("Cluster_E",LayerCluster_E) ;
    mytree->SetBranchAddress("Cluster_eta",LayerCluster_eta) ;
    mytree->SetBranchAddress("Cluster_phi",LayerCluster_phi) ;
    mytree->SetBranchAddress("Cluster_layer",LayerCluster_layer) ;
    mytree->SetBranchAddress("Cluster_mask_EM",LayerCluster_mask_EM) ;
    mytree->SetBranchAddress("Cluster_mask_HAD",LayerCluster_mask_HAD) ;
    mytree->SetBranchAddress("Cluster_mask_Trk",LayerCluster_mask_Trk) ;
    mytree->SetBranchAddress("Cluster_mask_MIP",LayerCluster_mask_MIP) ;
    //Branch for Particle
    TBranch *PFConstituents;
    std::vector<Double_t> *PFConstituents_eta;
    std::vector<Double_t> *PFConstituents_phi;
    std::vector<Double_t> *PFConstituents_E;
    std::vector<Int_t> *PFConstituents_charge;
    std::vector<Int_t> *PFConstituents_pdgID;
    mytree->SetBranchAddress("PFConstituents_eta", &PFConstituents_eta,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_phi", &PFConstituents_phi,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_charge", &PFConstituents_charge,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_E", &PFConstituents_E,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_pdgId",&PFConstituents_pdgID,&PFConstituents);
    
    Int_t nentries = (Int_t)mytree->GetEntries();
    //Setting for historgram 
    TH1D *Quark_barycentres = new TH1D("Barycentres for q"," Quark Barycentre" , 50 , 0, max_z-min_z) ;
    TH1D *Gluon_barycentres = new TH1D("Barycentres for g"," Gluon Barycentre" , 50 , 0, max_z-min_z) ;
    Quark_barycentres->SetFillColor(kRed); 
    Quark_barycentres->GetYaxis()->SetTitle("Events");
    Quark_barycentres->GetXaxis()->SetTitle("Depth");
    Quark_barycentres->GetXaxis()->SetTitleSize(font_size);
    Quark_barycentres->GetYaxis()->SetTitleSize(font_size);
    Quark_barycentres->GetXaxis()->SetTitleOffset(font_offset);
    Quark_barycentres->GetYaxis()->SetTitleOffset(font_offset);
    Quark_barycentres->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_barycentres->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_barycentres->SetLineColor(kBlack);

    Gluon_barycentres->SetFillColor(kGreen);
    Gluon_barycentres->GetYaxis()->SetTitle("Events");
    Gluon_barycentres->GetXaxis()->SetTitle("Depth");
    Gluon_barycentres->GetXaxis()->SetTitleSize(font_size);
    Gluon_barycentres->GetYaxis()->SetTitleSize(font_size);
    Gluon_barycentres->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_barycentres->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_barycentres->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_barycentres->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_barycentres->SetLineColor(kBlack);
    for( Int_t ev = 0 ; ev < nentries ; ev++)
    {
          mytree->GetEntry(ev);
          std::vector<Bool_t> LayerCluster_tagger;
          for(int j = 0 ; j < LayerCluster_track;j++)
          {
               LayerCluster_tagger.push_back(false);
          }
          for(int i=0; i<genparticle_track  ; i++)
          {
               bool quark=false;
               bool gluon=false;
               if(genparticle_pdgId[i]==Gluon_ID)
               {
                   gluon=true;
               }
               else if ((TMath::Abs(genparticle_pdgId[i])<=Quark_upper_ID&&TMath::Abs(genparticle_pdgId[i])>=Quark_lower_ID))
               {
                   quark=true; 
               }
               if((gluon||quark)&&genparticle_status[i] == 23)
               {
                   Double_t dR[ak4pfjets_track];
                   Double_t minimum_DeltaR=0;
                   int min_pos =0;
                   bool flag = false ;
                   for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                   {
                       dR[j] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[j] - genparticle_eta[i]),2) + TMath::Power((ak4pfjets_phi[j] - genparticle_phi[i]), 2 ) );
                   }
                   
                   Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos , flag );
                   if( minimum_DeltaR < 0.1 && TMath::Abs(ak4pfjets_eta[min_pos]) < hgcal_upper_eta \
                           && TMath::Abs( ak4pfjets_eta[min_pos])>hgcal_lower_eta  \
                           && ak4pfjets_pt[min_pos]>genparticle_pt[i]*constraint_coef_pt\
                           )
                   {
                       //Genparticle is Gluon or Quark
                       //Random choose one quark jet or gluon jet from events
                       Int_t particle_index;
                       if(min_pos >0 )
                       {
                           particle_index = ak4pfjets_tagging->at(min_pos-1);
                       }
                       else
                       {
                           particle_index=0;
                       }
                       Int_t particle_index_end = ak4pfjets_tagging->at(min_pos);
                   //Loop for particles 
                       for(Int_t k = 0 ; k < LayerCluster_track; k++ )
                       {
                           Double_t dR_=500;
                           int min_index=-1;
                           double deltaR=0;
                           for(int layer=0;layer<hgcallayer_max;layer++)
                           {
                               if(layer+1==LayerCluster_layer[k])
                               { 
                                   bool electron=false;
                                   bool photon=false;
                                   bool charged_Had=false;
                                   bool neutral_Had=false;
                                   if(LayerCluster_mask_EM[k]==0 && LayerCluster_mask_Trk[k]==0 && LayerCluster_mask_HAD[k]!=0\
                                           &&TMath::Abs(PFConstituents_pdgID->at(particle_index)) == 11)
                                   {
                                       electron=true;
                                   }
                                   else if(LayerCluster_mask_EM[k]==0 && LayerCluster_mask_Trk[k]!=0 && LayerCluster_mask_HAD[k]!=0\
                                           &&TMath::Abs(PFConstituents_pdgID->at(particle_index)) == 22)
                                   {
                                       photon=true;
                                   }
                                   else if(LayerCluster_mask_EM[k]==0 && LayerCluster_mask_Trk[k]==0 && LayerCluster_mask_HAD[k]==0\
                                           && PFConstituents_charge != 0
                                           )
                                   {
                                       charged_Had=true;
                                   }
                                   else if(LayerCluster_mask_EM[k]>0 && LayerCluster_mask_Trk[k]!=0 && LayerCluster_mask_HAD[k]==0\
                                           && PFConstituents_charge == 0)
                                   {
                                       neutral_Had=true;
                                   }
                                   if(electron || photon || charged_Had || neutral_Had)
                                   {
                                       deltaR = TMath::Power(ak4pfjets_eta[min_pos]-LayerCluster_eta[k],2)\
                                                +TMath::Power(ak4pfjets_phi[min_pos]-LayerCluster_phi[k],2);
                                       deltaR = TMath::Sqrt(deltaR);
                                       if(deltaR < 0.8)
                                       {
                                           LayerCluster_tagger.at(k)=true;
                                       }
                                   }

                               }
                           }
                       }

                       //Find out the total energy for each layer and the corresponding pattern in each layer
                       /*
                       double barycenter[hgcallayer_max]={0};
                       double Energy2[hgcallayer_max]={0};
                       double Depth[hgcallayer_max]={0};
                       */
                       double Barycenter=0.;
                       double Total_E2=0.;
                       for(int layer=0;layer<hgcallayer_max;layer++)
                       {
                           for(int  k = 0 ; k < LayerCluster_track ; k++)
                           {
                               if(LayerCluster_layer[k]==layer+1 && LayerCluster_tagger.at(k))
                               {
                                   double energy = LayerCluster_E[k];
                                   double depth = LayerCluster_z[k]-min_z;
                                   /*
                                   Energy2[layer] += energy;
                                   Depth[layer] = depth;    
                                   */
                                   Barycenter += energy*energy*depth;
                                   Total_E2 += energy*energy;
                                   LayerCluster_tagger.at(k)=false;
                               }
                           }
                       }
                       /*
                       for(int layer=0;layer<hgcallayer_max;layer++)
                       {
                           Energy2[layer] = Energy2[layer]*Energy2[layer];
                           Barycenter +=Energy2[layer]*Depth[layer];
                           Total_E2+=Energy2[layer];
                       }
                       */
                       Barycenter = Barycenter/Total_E2;
                       if(quark)
                           Quark_barycentres->Fill(Barycenter);
                       else if(gluon)
                           Gluon_barycentres->Fill(Barycenter);
                  }
               }
          }
   }

    c1->cd();
    Quark_barycentres->Draw();
    c1->SaveAs("Barycenre_Quark.pdf");
    c1->Clear();
    c2->cd();
    Gluon_barycentres->Draw();
    c2->SaveAs("Barycenre_Gluon.pdf");
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
