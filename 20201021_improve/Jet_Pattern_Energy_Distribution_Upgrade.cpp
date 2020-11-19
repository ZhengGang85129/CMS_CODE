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

const double canvas_width = 1200;
const double font_size = 0.055;
const double font_offset= 0.6;

const float Left_Margin = 0.12;
const float Right_Margin = 0.1;
const float Bottom_Margin = 0.12;
const float Top_Margin = 0.1;

const float marker_size = 0.8;
const int hgcallayer_max=50;
void Find_min(Double_t [] , Int_t , Int_t , Double_t & , int & ,bool &);
void Normalization(Double_t [] , const Int_t);

void Jet_Pattern_Energy_Distribution_Upgrade()
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
    
    Int_t hgcalclst_track;
    Int_t hgcalclst_layer[kMax];
    Int_t hgcalclst_size[kMax];
    Double_t hgcalclst_x[kMax];
    Double_t hgcalclst_y[kMax];
    Double_t hgcalclst_z[kMax];
    Double_t hgcalclst_E[kMax];
    Double_t hgcalclst_eta[kMax];
    Double_t hgcalclst_phi[kMax];

    mytree->SetBranchAddress("hgcalclst_track", &hgcalclst_track ) ;
    mytree->SetBranchAddress("hgcalclst_E",hgcalclst_E) ;
    mytree->SetBranchAddress("hgcalclst_eta",hgcalclst_eta) ;
    mytree->SetBranchAddress("hgcalclst_phi",hgcalclst_phi) ;
    mytree->SetBranchAddress("hgcalclst_layer",hgcalclst_layer) ;
    //Branch for Particle
    TBranch *PFConstituents;
    std::vector<Double_t> *PFConstituents_eta;
    std::vector<Double_t> *PFConstituents_phi;
    std::vector<Double_t> *PFConstituents_E;
    std::vector<Int_t> *PFConstituents_charge;
    mytree->SetBranchAddress("PFConstituents_eta", &PFConstituents_eta,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_phi", &PFConstituents_phi,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_charge", &PFConstituents_charge,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_E", &PFConstituents_E,&PFConstituents);

    /*Graph*/
    
    Int_t nentries = (Int_t)mytree->GetEntries();
    Double_t Energy_q[hgcallayer_max]={0};
    Double_t Energy_g[hgcallayer_max]={0};
    TH2D *Pattern_Quark[hgcallayer_max];
    TH2D *Pattern_Gluon[hgcallayer_max];
    TH2D *Jet_q = new TH2D("Quark Jet","Jet",25,0.,2.5,20,-3.,-1.);
    TH2D *Jet_g = new TH2D("Gluon Jet","Jet",25,-2.5,0.,20,-3.,-1);
    char name_q[50],name_g[50],title_q[50],title_g[50];
    sprintf(name_q,"Pattern_Quark");
    sprintf(name_g,"Pattern_Gluon");
    for(int index = 0 ; index < hgcallayer_max ;index++)
    {
        sprintf(name_q,"Pattern_Quark%d",index);
        sprintf(name_g,"Pattern_Gluon%d",index);
        sprintf(title_q,"Quark Jet Pattern in Layer%d",index+1);
        sprintf(title_g,"Gluon Jet Pattern in Layer%d",index+1);
        Pattern_Quark[index]=new TH2D(name_q,title_q,25,0.,2.5,20,-3.,-1.);
        Pattern_Gluon[index]=new TH2D(name_g,title_g,25,-2.5,0.,20,-3.,-1);  
        
        Pattern_Quark[index]->GetYaxis()->SetTitle("#eta");
        Pattern_Quark[index]->GetXaxis()->SetTitle("#phi");
        Pattern_Quark[index]->GetXaxis()->SetTitleSize(0.8*font_size);
        Pattern_Quark[index]->GetYaxis()->SetTitleSize(0.8*font_size); 
        Pattern_Quark[index]->GetXaxis()->SetTitleOffset(font_offset);
        Pattern_Quark[index]->GetYaxis()->SetTitleOffset(font_offset);

        Pattern_Gluon[index]->GetYaxis()->SetTitle("#eta");
        Pattern_Gluon[index]->GetXaxis()->SetTitle("#phi");
        Pattern_Gluon[index]->GetXaxis()->SetTitleSize(0.8*font_size);
        Pattern_Gluon[index]->GetYaxis()->SetTitleSize(0.8*font_size); 
        Pattern_Gluon[index]->GetXaxis()->SetTitleOffset(font_offset);
        Pattern_Gluon[index]->GetYaxis()->SetTitleOffset(font_offset);
    }
    Jet_q->GetYaxis()->SetTitle("#eta");
    Jet_q->GetXaxis()->SetTitle("#phi");
    Jet_q->GetXaxis()->SetTitleSize(0.8*font_size);
    Jet_q->GetYaxis()->SetTitleSize(0.8*font_size); 
    Jet_q->GetXaxis()->SetTitleOffset(font_offset);
    Jet_q->GetYaxis()->SetTitleOffset(font_offset);
    Jet_q->SetMarkerColor(kRed);
    Jet_q->SetMarkerSize(marker_size);
    Jet_q->SetMarkerStyle(3);

    Jet_g->GetYaxis()->SetTitle("#eta");
    Jet_g->GetXaxis()->SetTitle("#phi");
    Jet_g->GetXaxis()->SetTitleSize(0.8*font_size);
    Jet_g->GetYaxis()->SetTitleSize(0.8*font_size); 
    Jet_g->GetXaxis()->SetTitleOffset(font_offset);
    Jet_g->GetYaxis()->SetTitleOffset(font_offset);
    Jet_g->SetMarkerColor(kGreen);
    Jet_g->SetMarkerSize(marker_size);
    Jet_g->SetMarkerStyle(3);
    bool quark_one =true;
    bool gluon_one =true;
    int ctr_quark =0;
    int ctr_gluon =0;
    TGraph *Quark_EnergyDistribution = new TGraph() ;
    TGraph *Gluon_EnergyDistribution = new TGraph() ;
    for( Int_t ev = 0 ; ev < nentries ; ev++)
    {
          mytree->GetEntry(ev);
          std::vector<Bool_t> hgcalclst_tagger;
          for(int j = 0 ; j < hgcalclst_track;j++)
          {
               hgcalclst_tagger.push_back(false);
          }
          for(int i=0; i<genparticle_track  ; i++)
          {
               if(((TMath::Abs(genparticle_pdgId[i])<4&&TMath::Abs(genparticle_pdgId[i])>0)||genparticle_pdgId[i]==21)&&(genparticle_status[i] == 23))
               {
                   Double_t dR[ak4pfjets_track];
                   Double_t minimum_DeltaR=0;
                   int min_pos =0;
                   bool flag = false ;
                   bool quark=false;
                   bool gluon=false;
                   for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                   {
                       dR[j] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[j] - genparticle_eta[i]),2) + TMath::Power((ak4pfjets_phi[j] - genparticle_phi[i]), 2 ) );
                   }
                   
                   Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos , flag );
                   if( minimum_DeltaR < 0.1 && TMath::Abs(ak4pfjets_eta[min_pos]) < hgcal_upper_eta \
                           && TMath::Abs( ak4pfjets_eta[min_pos])>hgcal_lower_eta  \
                           && ak4pfjets_pt[min_pos]>genparticle_pt[i]*constraint_coef_pt)
                   {
                       if(genparticle_pdgId[i]==21)
                       {
                           gluon=true;
                           ctr_gluon++;
                       }
                       if((TMath::Abs(genparticle_pdgId[i])<4&&TMath::Abs(genparticle_pdgId[i])>0))
                       {
                           quark=true;
                           ctr_quark++;
                       }
                       if((ctr_quark==selected_q&&quark_one)||(ctr_gluon=selected_g&&gluon_one))
                       {
                           Int_t mindex;
                           if(gluon)
                           {
                               Jet_g->Fill(ak4pfjets_phi[min_pos],ak4pfjets_eta[min_pos]);
                               gluon_one=false;
                           }
                           else if(quark)
                           {
                               Jet_q->Fill(ak4pfjets_phi[min_pos],ak4pfjets_eta[min_pos]);
                               quark_one=false;
                           }
                           if(min_pos >0 )
                           {
                               mindex = ak4pfjets_tagging->at(min_pos-1);
                           }
                           else
                           {
                               mindex=0;
                           }
                           Int_t mindex_end = ak4pfjets_tagging->at(min_pos);
                           for(;mindex < mindex_end;mindex++)
                           {
                               for(Int_t k = 0 ; k < hgcalclst_track; k++ )
                               {
                                   Double_t dR_=500;
                                   int min_index=-1;
                                   double deltaR=0;
                                   for(int layer=0;layer<hgcallayer_max;layer++)
                                   {
                                       if(layer+1==hgcalclst_layer[k])
                                       {
                                           deltaR = TMath::Power(PFConstituents_eta->at(mindex)-hgcalclst_eta[k],2)+TMath::Power(PFConstituents_phi->at(mindex)-hgcalclst_phi[k],2);
                                           deltaR = TMath::Sqrt(deltaR);
                                           if(deltaR < dR_ && deltaR < part_dR_cut)
                                           {
                                               if(min_index>-1)
                                                   hgcalclst_tagger.at(min_index)=false;
                                               hgcalclst_tagger.at(k)=true;
                                               dR_=deltaR;
                                               min_index=k;
                                           }
                                       }
                                   }
                               }

                           }
                           for(int layer=0;layer<hgcallayer_max;layer++)
                           {
                               for(int  k = 0 ; k < hgcalclst_track ; k++)
                               {
                                   if(hgcalclst_layer[k]==layer+1 && hgcalclst_tagger.at(k))
                                   {
                                       if(quark)
                                       {
                                           Energy_q[layer] += hgcalclst_E[k];
                                           Pattern_Quark[layer]->Fill(hgcalclst_phi[k],hgcalclst_eta[k]);
                                           Int_t bin = Pattern_Quark[layer]->FindBin(hgcalclst_phi[k],hgcalclst_eta[k]);
                                           Pattern_Quark[layer]->SetBinContent(bin,hgcalclst_E[k]);
                                       }
                                       else if(gluon)
                                       {
                                           Energy_g[layer] += hgcalclst_E[k];
                                           Pattern_Gluon[layer]->Fill(hgcalclst_phi[k],hgcalclst_eta[k]);
                                           Int_t bin = Pattern_Gluon[layer]->FindBin(hgcalclst_phi[k],hgcalclst_eta[k]);
                                           Pattern_Gluon[layer]->SetBinContent(bin,hgcalclst_E[k]);
                                       }
                                       hgcalclst_tagger.at(k)=false;
                                   }
                               }
                           }
                       }
                  }
               }
          }
   }

/*Setting of Graph*/
    Quark_EnergyDistribution->Set(hgcallayer_max);
    Quark_EnergyDistribution->SetNameTitle("E Distribution for quark","Energy Distribution in Layers");
    Quark_EnergyDistribution->GetXaxis()->SetTitle("Layer number");
    Quark_EnergyDistribution->GetYaxis()->SetTitle("Total Energy per Layer");
    Quark_EnergyDistribution->GetXaxis()->SetTitleSize(0.8*font_size);
    Quark_EnergyDistribution->GetYaxis()->SetTitleSize(0.8*font_size); 
    Quark_EnergyDistribution->GetXaxis()->SetTitleOffset(font_offset);
    Quark_EnergyDistribution->GetYaxis()->SetTitleOffset(font_offset);
    Quark_EnergyDistribution->SetMarkerColor(kRed);
    
    Gluon_EnergyDistribution->Set(hgcallayer_max);
    Gluon_EnergyDistribution->SetNameTitle("E Distribution for Gluon","Energy Distribution in Layers");
    Gluon_EnergyDistribution->GetXaxis()->SetTitle("Layer number");
    Gluon_EnergyDistribution->GetYaxis()->SetTitle("Total Energy per Layer");
    Gluon_EnergyDistribution->GetXaxis()->SetTitleSize(font_size);
    Gluon_EnergyDistribution->GetYaxis()->SetTitleSize(font_size);
    Gluon_EnergyDistribution->GetYaxis()->SetTitleOffset(0.1*font_offset);
    Gluon_EnergyDistribution->GetXaxis()->SetTitleOffset(0.1*font_offset);
    Gluon_EnergyDistribution->SetMarkerColor(kGreen);  
    
    for(int i = 0; i < hgcallayer_max; i++)
    {
        Quark_EnergyDistribution->SetPoint(i+1,i+1,Energy_q[i]);
        Gluon_EnergyDistribution->SetPoint(i+1,i+1,Energy_g[i]);
    }

    c1->cd();
    Quark_EnergyDistribution->Draw("A*");
    c1->SaveAs("HGCal_Distribution_quark.pdf");
    c1->Clear();
    c2->cd();
    Gluon_EnergyDistribution->Draw("A*");
    c2->SaveAs("HGCal_Distribution_gluon.pdf");
    c2->Clear();
    for(int i = 0;i <hgcallayer_max;i++)
    {
        c1->cd();
        Pattern_Quark[i]->Draw("COLZ");
        c1->Modified();
        c1->Update();
        Jet_q->Draw("same");
        
        c2->cd();
        Pattern_Gluon[i]->Draw("COLZ");
        c2->Modified();
        c2->Update();
        Jet_g->Draw("same");
        if(i == 0)
        {
            c1->SaveAs("HGCal_Quark_Pattern.pdf(");
            c2->SaveAs("HGCal_Gluon_Pattern.pdf(");
        }
        else if(i==hgcallayer_max-1)
        {
            c1->SaveAs("HGCal_Quark_Pattern.pdf)");
            c2->SaveAs("HGCal_Gluon_Pattern.pdf)");
        }
        else
        {
            c1->SaveAs("HGCal_Quark_Pattern.pdf");
            c2->SaveAs("HGCal_Gluon_Pattern.pdf");
        }
        c1->Clear();
        c2->Clear();
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
