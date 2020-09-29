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

const double hgcal_lower_eta = 1.5;
const double hgcal_upper_eta = 3.;

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
void Find_min(Double_t [] , Int_t , Int_t , Double_t & , int & ,bool &);
void Normalization(Double_t [] , const Int_t);

void Matching_Hit_Jet()
{
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
    
    mytree->SetBranchAddress("ak4pfjets_track", &ak4pfjets_track );
    mytree->SetBranchAddress("ak4pfjets_charge", &ak4pfjets_charge[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_pdgId", &ak4pfjets_pdgId );
    mytree->SetBranchAddress("ak4pfjets_px", &ak4pfjets_px[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_py", &ak4pfjets_py[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_pz", &ak4pfjets_pz[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_pt", &ak4pfjets_pt[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_mass", &ak4pfjets_mass[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_mt", &ak4pfjets_mt[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_phi", &ak4pfjets_phi[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_theta", &ak4pfjets_theta[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_eta", &ak4pfjets_eta[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_rapidity", &ak4pfjets_rapidity[ak4pfjets_track] );
    mytree->SetBranchAddress("ak4pfjets_E", &ak4pfjets_E[ak4pfjets_track] );

    //Branch for HGCal 
    
    Int_t hgcalclst_track;
    Int_t hgcalclst_size[kMax];
    Double_t hgcalclst_x[kMax];
    Double_t hgcalclst_y[kMax];
    Double_t hgcalclst_z[kMax];
    Double_t hgcalclst_E[kMax];
    Double_t hgcalclst_eta[kMax];
    Double_t hgcalclst_phi[kMax];
    
    mytree->SetBranchAddress("hgcalclst_track", &hgcalclst_track ) ;
    mytree->SetBranchAddress("hgcalclst_x",hgcalclst_x) ;
    mytree->SetBranchAddress("hgcalclst_y",hgcalclst_y) ;
    mytree->SetBranchAddress("hgcalclst_z",hgcalclst_z) ;
    mytree->SetBranchAddress("hgcalclst_E",hgcalclst_E) ;
    mytree->SetBranchAddress("hgcalclst_eta",hgcalclst_eta) ;
    mytree->SetBranchAddress("hgcalclst_phi",hgcalclst_phi) ;

    TH2D *Quark_HGCalhits = new TH2D("HGCal Hits for q"," Quark hits " , 200 , -4, 4 , 200 , -4 , 4 ) ;
    TH2D *Gluon_HGCalhits = new TH2D("HGCal Hits for g"," Gluon hits " , 200 , -4, 4 , 200 , -4 , 4 ) ;
    TH2D *Quark_PFJets = new TH2D("PFJets for q ", " pfjets for q " , 200 , -4 , 4 , 200 , -4 , 4 );
    TH2D *Gluon_PFJets = new TH2D("PFJets for g ", " pfjets for g " , 200 , -4 , 4 , 200 , -4 , 4 );
    TLine *boundary1 =new TLine( -4. , 1.5 , 4. , 1.5 ); 
    TLine *boundary2 =new TLine( -4. , 3. , 4. , 3. );
    TLine *boundary3 =new TLine( -4. , -1.5 , 4. , -1.5 ); 
    TLine *boundary4 =new TLine( -4. , -3. , 4. , -3. );

/*Setting of objects*/
    Quark_HGCalhits->GetXaxis()->SetTitle("#phi");
    Quark_HGCalhits->GetYaxis()->SetTitle("#eta");
    Quark_HGCalhits->GetXaxis()->SetTitleSize(font_size);
    Quark_HGCalhits->GetYaxis()->SetTitleSize(font_size); 
    Quark_HGCalhits->GetXaxis()->SetTitleOffset(font_offset);
    Quark_HGCalhits->GetYaxis()->SetTitleOffset(font_offset);

    Quark_PFJets->SetMarkerColor(kRed);
    Quark_PFJets->SetMarkerSize(marker_size);
    Quark_PFJets->SetMarkerStyle(3);
    
    Gluon_HGCalhits->GetXaxis()->SetTitle("#phi");
    Gluon_HGCalhits->GetYaxis()->SetTitle("#eta");
    Gluon_HGCalhits->GetXaxis()->SetTitleSize(font_size);
    Gluon_HGCalhits->GetYaxis()->SetTitleSize(font_size);
    Gluon_HGCalhits->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_HGCalhits->GetXaxis()->SetTitleOffset(font_offset);
    
    Gluon_PFJets->SetMarkerColor(kGreen);
    Gluon_PFJets->SetMarkerSize(marker_size);
    Gluon_PFJets->SetMarkerStyle(3); 
    
    boundary1->SetLineColor(kBlue);
    boundary2->SetLineColor(kBlue);
    boundary3->SetLineColor(kBlue); 
    boundary4->SetLineColor(kBlue);
/*Setting of Histogram*/
    
    Int_t nentries = (Int_t)mytree->GetEntries();
    int quark_ctr = 0 ;
    int gluon_ctr = 0 ;
    for( Int_t ev = 0 ; ev < nentries ; ev++)
    {
       mytree->GetEntry(ev);
      for(int i = 0; i < genparticle_track  ; i ++)
      {
       if(TMath::Abs(genparticle_pdgId[ i ])< 4 && TMath::Abs(genparticle_pdgId[ i ]) > 0 && genparticle_status[ i ] == 23  && quark_ctr < Number_of_HIST)
       {
           Double_t dR[ ak4pfjets_track ] ;
           Double_t minimum_DeltaR;
           int min_pos ;
           bool flag = false ;
           for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                  {
                      dR[ j ] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[ j ] - genparticle_eta[ i ]),2) + TMath::Power((ak4pfjets_phi[ j ] - genparticle_phi[ i ]), 2 ) );
                    
                  }
           
           Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos , flag );
      //     Quark_min_DeltaR->Fill( minimum_DeltaR );
           if( minimum_DeltaR < 0.1 && TMath::Abs(ak4pfjets_eta[ min_pos ]) < hgcal_upper_eta \
                   && TMath::Abs( ak4pfjets_eta[ min_pos ]) >  hgcal_lower_eta  \
                   && ak4pfjets_pt[  min_pos ] > genparticle_pt[ min_pos ] * constraint_coef_pt )
           {
               Double_t dR_; 
               quark_ctr++;
               Quark_PFJets->Fill(ak4pfjets_phi[ min_pos ],ak4pfjets_eta[ min_pos ]);
               for(Int_t k = 0 ; k < hgcalclst_track ; k ++ )
                    {
                        dR_ = TMath::Sqrt(TMath::Power(hgcalclst_eta[ k ] - ak4pfjets_eta[ min_pos ] , 2 ) + TMath::Power( hgcalclst_phi[ k ] - ak4pfjets_phi[ min_pos ] , 2 ) );
                        if( dR_ < R_cut )
                        {
                            Quark_HGCalhits->Fill(hgcalclst_phi[k],hgcalclst_eta[k]);
                            Int_t bin = Quark_HGCalhits->FindBin(hgcalclst_phi[k],hgcalclst_eta[k]);
                            Quark_HGCalhits->SetBinContent(bin,hgcalclst_E[k]);
                            //printf("%f\n",Quark_HGCalhits->GetXaxis()->GetBinWidth(bin));
                            //printf("%f\n",Quark_HGCalhits->GetYaxis()->GetBinWidth(bin));
                            //HGCal_ZX_profile_q->Fill( hgcalclst_z[ k ] , hgcalclst_x[ k ]);
                        }

                    }
               c1->cd();
               Quark_HGCalhits->GetZaxis()->SetTickSize(0.01);
               Quark_HGCalhits->Draw("COLZ") ; 
               c1->Modified();
               c1->Update();
               boundary1->Draw("same");
               boundary2->Draw("same");
               boundary3->Draw("same");
               boundary4->Draw("same");
               Quark_PFJets->Draw("same");
               if(quark_ctr == 1 )
                   c1->SaveAs("HGCal_Distribution_quark.pdf(");
               else
                   c1->SaveAs("HGCal_Distribution_quark.pdf");
               Quark_HGCalhits->Reset("ICESM");
               Quark_PFJets->Reset("ICESM");
               c1->Clear();
           }
       }
       else if(genparticle_pdgId[ i ] == 21 && genparticle_status[ i ] == 23 && gluon_ctr < Number_of_HIST)
       {
           Double_t dR[ ak4pfjets_track ] ;
           Double_t minimum_DeltaR ;
           bool flag = false ; 
           int min_pos ;
           for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                  {
                      dR[ j ] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[ j ] - genparticle_eta[ i ]),2) + TMath::Power((ak4pfjets_phi[ j ] - genparticle_phi[ i ]), 2 ) );
                  }
           Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos ,flag);
           if( minimum_DeltaR < 0.1  && TMath::Abs(ak4pfjets_eta[ min_pos ]) < hgcal_upper_eta && TMath::Abs( ak4pfjets_eta[ min_pos ]) >  hgcal_lower_eta && ak4pfjets_pt[  min_pos ] > \
                   genparticle_pt[ min_pos ] * constraint_coef_pt)
           {
               Double_t dR_ ;
               Gluon_PFJets->Fill(ak4pfjets_phi[ min_pos ],ak4pfjets_eta[ min_pos]);
               gluon_ctr++;
               for(Int_t k = 0 ; k < hgcalclst_track ; k ++)
               {
                   dR_ = TMath::Sqrt(TMath::Power(hgcalclst_eta[ k ] - ak4pfjets_eta[ min_pos ] , 2 ) + TMath::Power( hgcalclst_phi[ k ] - ak4pfjets_phi[ min_pos ] , 2 ) );
                   if( dR_ < R_cut)
                   {
                       Gluon_HGCalhits->Fill( hgcalclst_phi[ k ] , hgcalclst_eta[ k ] );
                       Int_t bin = Gluon_HGCalhits->FindBin(hgcalclst_phi[k],hgcalclst_eta[k]);
                       Gluon_HGCalhits->SetBinContent(bin,hgcalclst_E[k]);
//                       HGCal_ZX_profile_g->Fill( hgcalclst_z[ k ] , hgcalclst_x[ k ]);  
                   }
               }
               c2->cd();
               c2->Modified();
               Gluon_HGCalhits->Draw("COLZ") ;
               c2->Update();
               boundary1->Draw("same");
               boundary2->Draw("same");
               boundary3->Draw("same");
               boundary4->Draw("same");
               Gluon_PFJets->Draw("same"); 
               if(quark_ctr == 1 )
                   c2->SaveAs("HGCal_Distribution_gluon.pdf(");
               else
                   c2->SaveAs("HGCal_Distribution_gluon.pdf");
               Gluon_PFJets->Reset("ICESM");  
               Gluon_HGCalhits->Reset("ICESM");
               c2->Clear();
           }
       }
      }
   }
    c1->SaveAs("HGCal_Distribution_quark.pdf)");
    c2->SaveAs("HGCal_Distribution_gluon.pdf)");
    c1->Clear();
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
void Normalization(Double_t arr[] , const Int_t size )
{
    Double_t avg = 0 ; 
    Double_t std = 0;
    for(int i=0 ; i < size ; i++)
    {
        avg+=arr[i];
        std+=arr[i]*arr[i];
    }
    avg /=static_cast< double >(size);
    std /=static_cast< double >(size);
    std -= avg;
    std = TMath::Sqrt(std);
    for(int i=0 ; i < size ; i++)
    {
        arr[i]=(arr[i]-avg)/std;
    }

}
