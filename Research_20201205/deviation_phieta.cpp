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
const Double_t constraint_coef_pt = .8;

const int Number_of_HIST= 10;
const double step_of_value=0.3;
const double Init_value=0.3;

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

void deviation_phieta()
{
    delete gROOT->FindObject("Radiance_31");
    TSeqCollection* canvases = gROOT->GetListOfCanvases();
    TIter next(gROOT->GetListOfCanvases());
    while(c = (TCanvas*)next())
    {
        delete c;
    }
    gStyle->SetOptStat(0);
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
    TH1D *Quark_phi_deviation[Number_of_HIST];
    TH1D *Gluon_phi_deviation[Number_of_HIST];
    TH1D *Quark_eta_deviation[Number_of_HIST];
    TH1D *Gluon_eta_deviation[Number_of_HIST];
    TH1D *Quark_dr_deviation[Number_of_HIST];
    TH1D *Gluon_dr_deviation[Number_of_HIST];
    TLegend legend_phi[Number_of_HIST];
    TLegend legend_eta[Number_of_HIST];
    TLegend legend_dr[Number_of_HIST];
    char titleg_phi[80];
    char titleq_phi[80];
    char titleq_eta[80];
    char titleg_eta[80];
    char titleq_dr[80];
    char titleg_dr[80];

    char nameq_phi[80];
    char nameg_phi[80];
    char nameq_eta[80];
    char nameg_eta[80];
    char nameg_dr[80];
    char nameq_dr[80];
    
    for(int i = 0 ;i <Number_of_HIST; i++)
    {
        double para=0.5;
        sprintf(titleg_phi,"#it{#phi} deviation with E^{%.2f}",para+static_cast<double>(i)*0.5);
        sprintf(titleq_phi,"#it{#phi} deviation with E^{%.2f}",para+static_cast<double>(i)*0.5);
        sprintf(titleg_eta,"#eta deviation with E^{%.2f}",para+static_cast<double>(i)*0.5);
        sprintf(titleq_eta,"#eta deviation with E^{%.2f}",para+static_cast<double>(i)*0.5);
        sprintf(titleg_dr,"#DeltaR deviation with E^{%.2f}",para+static_cast<double>(i)*0.5);
        sprintf(titleq_dr,"#DeltaR with E^{%.2f}",para+static_cast<double>(i)*0.5);
        
        sprintf(nameq_phi,"#phi deviation with E^{%.2f} for q",para+static_cast<double>(i)*0.5);
        sprintf(nameg_phi,"#phi deviation with E^{%.2f} for g",para+static_cast<double>(i)*0.5);
        sprintf(nameq_eta,"#eta deviation with E^{%.2f} for q",para+static_cast<double>(i)*0.5);
        sprintf(nameg_eta,"#eta deviation with E^{%.2f} for g",para+static_cast<double>(i)*0.5);
        sprintf(nameq_dr,"#Delta R deviation with E^{%.2f} for q",para+static_cast<double>(i)*0.5);
        sprintf(nameg_dr,"#Delta R deviation with E^{%.2f} for g",para+static_cast<double>(i)*0.5);
        
        Quark_phi_deviation[i]=new TH1D(nameq_phi,titleq_phi,50,0.,0.3);
        Quark_phi_deviation[i]->SetFillColorAlpha(kRed,0.15);
        Quark_phi_deviation[i]->GetYaxis()->SetTitle("Events");
        Quark_phi_deviation[i]->GetXaxis()->SetTitle("Unit(#it{#phi})");
        Quark_phi_deviation[i]->GetXaxis()->SetTitleSize(font_size);
        Quark_phi_deviation[i]->GetYaxis()->SetTitleSize(font_size);
        Quark_phi_deviation[i]->GetXaxis()->SetTitleOffset(font_offset);
        Quark_phi_deviation[i]->GetYaxis()->SetTitleOffset(font_offset);
        Quark_phi_deviation[i]->GetXaxis()->SetLabelSize(0.4*font_size);
        Quark_phi_deviation[i]->GetYaxis()->SetLabelSize(0.4*font_size);
        Quark_phi_deviation[i]->SetLineWidth(2);
        Quark_phi_deviation[i]->SetLineColor(kRed);

        Gluon_phi_deviation[i]=new TH1D(nameg_phi,titleg_phi,50,0.,0.3);
        Gluon_phi_deviation[i]->SetFillColorAlpha(kBlue,0.15);
        Gluon_phi_deviation[i]->GetYaxis()->SetTitle("Events");
        Gluon_phi_deviation[i]->GetXaxis()->SetTitle("Unit(#it{#phi})");
        Gluon_phi_deviation[i]->GetXaxis()->SetTitleSize(font_size);
        Gluon_phi_deviation[i]->GetYaxis()->SetTitleSize(font_size);
        Gluon_phi_deviation[i]->GetXaxis()->SetTitleOffset(font_offset);
        Gluon_phi_deviation[i]->GetYaxis()->SetTitleOffset(font_offset);
        Gluon_phi_deviation[i]->GetXaxis()->SetLabelSize(0.4*font_size);
        Gluon_phi_deviation[i]->GetYaxis()->SetLabelSize(0.4*font_size);
        Gluon_phi_deviation[i]->SetLineWidth(2);
        Gluon_phi_deviation[i]->SetLineColor(kBlue);

        Quark_eta_deviation[i]=new TH1D(nameq_eta,titleq_eta,50,0.,0.3);
        Quark_eta_deviation[i]->SetFillColorAlpha(kRed,0.15);
        Quark_eta_deviation[i]->GetYaxis()->SetTitle("Events");
        Quark_eta_deviation[i]->GetXaxis()->SetTitle("Unit(#it{#eta})");
        Quark_eta_deviation[i]->GetXaxis()->SetTitleSize(font_size);
        Quark_eta_deviation[i]->GetYaxis()->SetTitleSize(font_size);
        Quark_eta_deviation[i]->GetXaxis()->SetTitleOffset(font_offset);
        Quark_eta_deviation[i]->GetYaxis()->SetTitleOffset(font_offset);
        Quark_eta_deviation[i]->GetXaxis()->SetLabelSize(0.4*font_size);
        Quark_eta_deviation[i]->GetYaxis()->SetLabelSize(0.4*font_size);
        Quark_eta_deviation[i]->SetLineWidth(2);
        Quark_eta_deviation[i]->SetLineColor(kRed);
        
        Gluon_eta_deviation[i]=new TH1D(nameg_eta,titleg_eta,50,0.,0.3);
        Gluon_eta_deviation[i]->SetFillColorAlpha(kBlue,0.15);
        Gluon_eta_deviation[i]->GetYaxis()->SetTitle("Events");
        Gluon_eta_deviation[i]->GetXaxis()->SetTitle("Unit(#it{#eta})");
        Gluon_eta_deviation[i]->GetXaxis()->SetTitleSize(font_size);
        Gluon_eta_deviation[i]->GetYaxis()->SetTitleSize(font_size);
        Gluon_eta_deviation[i]->GetXaxis()->SetTitleOffset(font_offset);
        Gluon_eta_deviation[i]->GetYaxis()->SetTitleOffset(font_offset);
        Gluon_eta_deviation[i]->GetXaxis()->SetLabelSize(0.4*font_size);
        Gluon_eta_deviation[i]->GetYaxis()->SetLabelSize(0.4*font_size);
        Gluon_eta_deviation[i]->SetLineWidth(2);
        Gluon_eta_deviation[i]->SetLineColor(kBlue);
        
        Quark_dr_deviation[i]=new TH1D(nameq_dr,titleq_dr,50,0.,0.30);
        Quark_dr_deviation[i]->SetFillColorAlpha(kRed,0.15);
        Quark_dr_deviation[i]->GetYaxis()->SetTitle("Events");
        Quark_dr_deviation[i]->GetXaxis()->SetTitle("Unit(#it{#DetlaR})");
        Quark_dr_deviation[i]->GetXaxis()->SetTitleSize(font_size);
        Quark_dr_deviation[i]->GetYaxis()->SetTitleSize(font_size);
        Quark_dr_deviation[i]->GetXaxis()->SetTitleOffset(font_offset);
        Quark_dr_deviation[i]->GetYaxis()->SetTitleOffset(font_offset);
        Quark_dr_deviation[i]->GetXaxis()->SetLabelSize(0.4*font_size);
        Quark_dr_deviation[i]->GetYaxis()->SetLabelSize(0.4*font_size);
        Quark_dr_deviation[i]->SetLineWidth(2);
        Quark_dr_deviation[i]->SetLineColor(kRed);
        
        Gluon_dr_deviation[i]=new TH1D(nameg_dr,titleg_dr,50,0.,0.30);
        Gluon_dr_deviation[i]->SetFillColorAlpha(kBlue,0.15);
        Gluon_dr_deviation[i]->GetYaxis()->SetTitle("Events");
        Gluon_dr_deviation[i]->GetXaxis()->SetTitle("Unit(#it{#DeltaR})");
        Gluon_dr_deviation[i]->GetXaxis()->SetTitleSize(font_size);
        Gluon_dr_deviation[i]->GetYaxis()->SetTitleSize(font_size);
        Gluon_dr_deviation[i]->GetXaxis()->SetTitleOffset(font_offset);
        Gluon_dr_deviation[i]->GetYaxis()->SetTitleOffset(font_offset);
        Gluon_dr_deviation[i]->GetXaxis()->SetLabelSize(0.4*font_size);
        Gluon_dr_deviation[i]->GetYaxis()->SetLabelSize(0.4*font_size);
        Gluon_dr_deviation[i]->SetLineWidth(2);
        Gluon_dr_deviation[i]->SetLineColor(kBlue);
       
        legend_phi[i].SetX1NDC(0.65);
        legend_phi[i].SetY1NDC(0.7);
        legend_phi[i].SetX2NDC(0.88);
        legend_phi[i].SetY2NDC(0.88);
        legend_phi[i].SetHeader("Representation");
        legend_phi[i].AddEntry(Quark_phi_deviation[i],"Quark #it{#phi} deviation","f");
        legend_phi[i].AddEntry(Gluon_phi_deviation[i],"Gluon #it{#phi} deviation","f");
        
        legend_eta[i].SetX1NDC(0.65);
        legend_eta[i].SetY1NDC(0.7);
        legend_eta[i].SetX2NDC(0.88);
        legend_eta[i].SetY2NDC(0.88);
        legend_eta[i].SetHeader("Representation");
        legend_eta[i].AddEntry(Quark_eta_deviation[i],"Quark #it{#eta} deviation","f");
        legend_eta[i].AddEntry(Gluon_eta_deviation[i],"Gluon #it{#eta} deviation","f");
        
        legend_dr[i].SetX1NDC(0.65);
        legend_dr[i].SetY1NDC(0.7);
        legend_dr[i].SetX2NDC(0.88);
        legend_dr[i].SetY2NDC(0.88);
        legend_dr[i].SetHeader("Representation");
        legend_dr[i].AddEntry(Quark_dr_deviation[i],"Quark #it{#DeltaR} deviation","f");
        legend_dr[i].AddEntry(Gluon_dr_deviation[i],"Gluon #it{#DeltaR} deviation","f");
    }
//Legend/
    
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
                       for(;particle_index<particle_index_end;particle_index++)
                       {
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
                       }

                       //Find out the total energy for each layer and the corresponding pattern in each layer
                       double Total_E_power[Number_of_HIST]={0.};
                       double dev_phi[Number_of_HIST] ={0.};//deviation of phi weigted by power of energy
                       double dev_eta[Number_of_HIST] = {0.};
                       double dev_dr[Number_of_HIST]={0.};
                       double avg_phi=0.;//Average phi
                       double avg_eta=0.;//Average eta
                       double avg_dr=0.;//Average delta r
                       int  n_ALC =0;//number of Associated Cluster
                       for(int layer=0;layer<hgcallayer_max;layer++)
                       {
                           for(int  k = 0 ; k < LayerCluster_track ; k++)
                           {
                               if(LayerCluster_layer[k]==layer+1 && LayerCluster_tagger.at(k))
                               {
                                   avg_phi+=LayerCluster_phi[k];
                                   avg_eta+=LayerCluster_eta[k];
                                   for(int number = 0 ; number < Number_of_HIST ; number++)
                                   {
                                       Total_E_power[number]+=TMath::Power(LayerCluster_E[k],Init_value+static_cast<double>(number)*step_of_value);
                                   }
                                   n_ALC++;
                               }
                           }
                       }
                       avg_phi/=static_cast<double>(n_ALC);
                       avg_eta/=static_cast<double>(n_ALC);
                       for(int layer=0;layer<hgcallayer_max;layer++)
                       {
                           for(int  k = 0 ; k < LayerCluster_track ; k++)
                           {
                               if(LayerCluster_layer[k]==layer+1 && LayerCluster_tagger.at(k))
                               {
                                   
                                   for(int number = 0 ; number < Number_of_HIST ; number++)
                                   {
                                       double phi_weight = TMath::Power(LayerCluster_E[k],static_cast<double>(number)*step_of_value+Init_value)*TMath::Abs(LayerCluster_phi[k]-avg_phi);
                                       double eta_weight = TMath::Power(LayerCluster_E[k],static_cast<double>(number)*step_of_value+Init_value)*TMath::Abs(LayerCluster_eta[k]-avg_eta); 
                                       dev_phi[number] += phi_weight;
                                       dev_eta[number] += eta_weight;
                                       dev_dr[number] += TMath::Sqrt((eta_weight*eta_weight+phi_weight*phi_weight));
                                   }
                                   //---Clean tagger---//
                                   LayerCluster_tagger.at(k)=false;
                               }
                           }
                       }

                       for(int number = 0 ; number < Number_of_HIST ; number++)
                       {
                           dev_phi[number]/=Total_E_power[number];
                           dev_eta[number]/=Total_E_power[number];
                           dev_dr[number]=TMath::Sqrt(dev_phi[number]*dev_phi[number]+dev_eta[number]*dev_eta[number]);
                           if(quark)
                           {
                               Quark_phi_deviation[number]->Fill(dev_phi[number]);
                               Quark_eta_deviation[number]->Fill(dev_eta[number]);
                               Quark_dr_deviation[number]->Fill(dev_dr[number]);
                           }
                           else if(gluon)
                           {
                               Gluon_phi_deviation[number]->Fill(dev_phi[number]);
                               Gluon_eta_deviation[number]->Fill(dev_eta[number]);
                               Gluon_dr_deviation[number]->Fill(dev_dr[number]);
                           }
                       }
                  }
               }
          }
   }
    char q_title_phi[80] ="";
    char q_title_eta[80] ="";
    char q_title_dr[80] = "";
    for(int number = 0 ; number<Number_of_HIST ; number++)
    {
        if(number==0)
        {
            sprintf(q_title_phi,"Phi_Deviation_pwr.pdf(");
            sprintf(q_title_eta,"Eta_Deviation_pwr.pdf(");
            sprintf(q_title_dr,"DeltaR_Deviation_pwr.pdf(");
        }
        else if(number==Number_of_HIST-1)
        {
            sprintf(q_title_phi,"Phi_Deviation_pwr.pdf)");
            sprintf(q_title_eta,"Eta_Deviation_pwr.pdf)");
            sprintf(q_title_dr,"DeltaR_Deviation_pwr.pdf)");

        }
        else 
        {
            sprintf(q_title_phi,"Phi_Deviation_pwr.pdf");
            sprintf(q_title_eta,"Eta_Deviation_pwr.pdf");
            sprintf(q_title_dr,"DeltaR_Deviation_pwr.pdf");

        }
        c1->cd();
        
        double scale_q = 1./Quark_phi_deviation[number]->Integral();
        double scale_g = 1./Gluon_phi_deviation[number]->Integral();
        
        Quark_phi_deviation[number]->Scale(scale_q);
        Gluon_phi_deviation[number]->Scale(scale_g);
        Gluon_phi_deviation[number]->SetMaximum(1.5 * Quark_phi_deviation[number]->GetMaximum());
        Gluon_phi_deviation[number]->Draw();
        Quark_phi_deviation[number]->Draw("SAME");
        legend_phi[number].Draw();
        
        c1->SaveAs(q_title_phi);
        c1->Clear();
        c2->cd(); 
        double scale_q1 = 1./Quark_eta_deviation[number]->Integral();
        double scale_g1 = 1./Gluon_eta_deviation[number]->Integral();
        
        Quark_eta_deviation[number]->Scale(scale_q1);
        Gluon_eta_deviation[number]->Scale(scale_g1);
        Gluon_eta_deviation[number]->SetMaximum(1.5 * Quark_eta_deviation[number]->GetMaximum());
        Gluon_eta_deviation[number]->Draw();
        Quark_eta_deviation[number]->Draw("SAME");
        legend_eta[number].Draw();
        
        c2->SaveAs(q_title_eta);
        c2->Clear();
        
        c1->cd();
        scale_q = 1./Quark_dr_deviation[number]->Integral();
        scale_g = 1./Gluon_dr_deviation[number]->Integral();
        Quark_dr_deviation[number]->Scale(scale_q);
        Gluon_dr_deviation[number]->Scale(scale_g);
        Gluon_dr_deviation[number]->SetMaximum(1.5 * Quark_dr_deviation[number]->GetMaximum());
        Gluon_dr_deviation[number]->Draw();
        Quark_dr_deviation[number]->Draw("SAME");
        legend_dr[number].Draw();
        
        c1->SaveAs(q_title_dr);
        c1->Clear();
    }
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
