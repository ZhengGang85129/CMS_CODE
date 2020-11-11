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
const Double_t pt_min = 30.;
const double canvas_width = 1200;
const double font_size = .035;
const double font_offset= .6;

const float Left_Margin = .12;
const float Right_Margin = .1;
const float Bottom_Margin = .12;
const float Top_Margin = .1;

const float marker_size = .8;
void Find_min(Double_t [] , Int_t , Int_t , Double_t & , int & ,bool &);
void Jet_Width()
{
    gROOT->Reset();
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
    const double max_sigma = 0.4;
    TH1D *Quark_width = new TH1D("Quark width","Average width for Quark" , 30 , 0, max_sigma ) ;

//For Gluon    
    TH1D *Gluon_width = new TH1D("Gluon width","Average width for Gluon" , 30 , 0, max_sigma) ;

/*Setting of objects*/

//For Quark   
    Quark_width->GetXaxis()->SetTitle("#sigma");
    Quark_width->GetYaxis()->SetTitle("Events");
    Quark_width->GetXaxis()->SetTitleSize(font_size);
    Quark_width->GetYaxis()->SetTitleSize(font_size); 
    Quark_width->GetXaxis()->SetTitleOffset(font_offset);
    Quark_width->GetYaxis()->SetTitleOffset(font_offset);
    Quark_width->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_width->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_width->SetFillColor(kRed);
    Quark_width->SetLineColor(kBlack);

//For Gluon
    Gluon_width->GetXaxis()->SetTitle("#sigma");
    Gluon_width->GetYaxis()->SetTitle("Events");
    Gluon_width->GetXaxis()->SetTitleSize(font_size);
    Gluon_width->GetYaxis()->SetTitleSize(font_size); 
    Gluon_width->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_width->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_width->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_width->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_width->SetFillColor(kGreen);
    Gluon_width->SetLineColor(kBlack);

    Int_t nentries = (Int_t)mytree->GetEntries();
    for( Int_t ev = 0 ; ev < nentries ; ev++)
    {
            mytree->GetEntry(ev);
            for(int i = 0; i < genparticle_track  ; i ++)
            {
                    if(genparticle_status[ i ]==23)
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
                                   && ak4pfjets_pt[min_pos] > genparticle_pt[ min_pos ]*constraint_coef_pt\
                                   && pt_min < ak4pfjets_pt[min_pos]
                                   )
                           {
                               double sq_pt = 0;
                               int part_total_number=0;
                               double M11 = 0;
                               double M12 = 0;
                               double M22 = 0;
                               Int_t mindex;
                               Int_t mindex_end = ak4pfjets_tagging->at(min_pos);
                               if(min_pos >0 )
                               {
                                   mindex = ak4pfjets_tagging->at(min_pos-1);
                               }
                               else
                               {
                                   mindex = 0;
                               }
                               for(;mindex < mindex_end;mindex++)
                               {
                                    if(PFConstituents_charge->at(mindex) !=0 ||(PFConstituents_charge->at(mindex) ==0\
                                                && PFConstituents_pt->at(mindex) >1))
                                    {
                                        sq_pt+=TMath::Power(PFConstituents_pt->at(mindex),2);
                                        part_total_number++;
                                        double delta_eta;
                                        double delta_phi;
                                        double pt_2;
                                        delta_eta = PFConstituents_eta->at(mindex)-ak4pfjets_eta[min_pos];
                                        delta_phi = PFConstituents_phi->at(mindex)-ak4pfjets_phi[min_pos];
                                        pt_2 = TMath::Power(PFConstituents_pt->at(mindex),2);
                                        M11 += TMath::Power(delta_eta,2) * pt_2;
                                        M22 += TMath::Power(delta_phi,2) * pt_2;
                                        M12 += -delta_eta*delta_phi*pt_2;
                                    }
                               }
                               Double_t width_1 , width_2  , width ;
                               width_1 = ((M11+M22)+TMath::Sqrt((M11-M22)*(M11-M22)+4.*M12*M12))/2.;
                               width_2 = ((M11-M22)+TMath::Sqrt((M11-M22)*(M11-M22)+4.*M12*M12))/2.;
                               width_1 /= sq_pt;
                               width_2 /= sq_pt;
                               width_1 = TMath::Sqrt(width_1);
                               width_2 = TMath::Sqrt(width_2);
                               width = TMath::Sqrt(width_1 * width_1 + width_2 * width_2 ); 
                               if(TMath::Abs(genparticle_pdgId[i])< 4 && TMath::Abs(genparticle_pdgId[i]) > 0)
                               {
                                   Quark_width->Fill(width);
                               }
                               else if(genparticle_pdgId[i]==21)
                               {
                                   Gluon_width->Fill(width);
                               }
                           }
                    }
                }
    }
    c1->cd();
    Quark_width->Draw();
    c1->SaveAs("width_quark.pdf");
    c1->Clear();
    c2->cd();
    Gluon_width->Draw();
    c2->SaveAs("width_gluon.pdf");
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

