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

using namespace std;
const Int_t kMax = 10000;

void read_genparticle();
void read_AK4PFJet();

void read()
{
    read_genparticle();
    read_AK4PFJet();
}

void read_genparticle()
{
//Branch for genparticle
    TFile *f = new TFile("/Users/chenzhenggang/Desktop/local_analysis/Ntuple.root");
    //TFile *f = new TFile("/Users/chenzhenggang/Desktop/local_analysis/Ntuple_1_event.root");
    TTree *mytree = (TTree*)f->Get("Tree");
    TCanvas *c1=new TCanvas("c1","c1",1600,1600);
    Int_t genparticle_track;
    Int_t genparticle_charge[kMax];
    Int_t genparticle_pdgId[kMax];
    Double32_t genparticle_px[kMax];
    Double32_t genparticle_py[kMax];
    Double32_t genparticle_pz[kMax];
    Double32_t genparticle_pt[kMax];
    Double32_t genparticle_mass[kMax];
    Double32_t genparticle_mt[kMax];
    Double32_t genparticle_phi[kMax];
    Double32_t genparticle_theta[kMax];
    Double32_t genparticle_eta[kMax];
    Double32_t genparticle_rapidity[kMax];
    
    mytree->SetBranchAddress("genparticle_track",&genparticle_track);
    mytree->SetBranchAddress("genparticle_charge",&genparticle_charge[genparticle_track]);
    mytree->SetBranchAddress("genparticle_pdgId",&genparticle_pdgId);
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
    
    TH1F *hgenparticle_track = new TH1F("gen_track"," track distribution ", 36, 0 , 3500);
    TH1F *hgenparticle_charge = new TH1F("gen_charge","charge distribution", 5, -2 , 2);
    TH1F *hgenparticle_pdgId = new TH1F("gen_pdgId","PDG_ID distribution", 11, 0 , 10);
    TH1F *hgenparticle_px = new TH1F("gen_px","Px distribution", 81, -4000, 4000 );
    TH1F *hgenparticle_py = new TH1F("gen_py","Py distribution", 161, -4000, 4000 );
    TH1F *hgenparticle_pz = new TH1F("gen_pz","Pz distribution", 130, -8000, 8000 );
    TH1F *hgenparticle_pt = new TH1F("gen_pt","Pt distribution", 21, 0, 2000 );
    TH1F *hgenparticle_mass = new TH1F("gen_mass","Mass distribution", 31, 0., 6. );
    TH1F *hgenparticle_mt = new TH1F("gen_mt","Transverse mass distribution", 41, 0, 4000 );
    TH1F *hgenparticle_phi = new TH1F("gen_phi","Phi distribution", 7, -TMath::Pi(),TMath::Pi() );
    TH1F *hgenparticle_theta = new TH1F("gen_theta","Theta distribution", 13, 0,  TMath::Pi());
    TH1F *hgenparticle_eta = new TH1F("gen_eta","Eta distribution", 61, -30000, 30000 );
    TH1F *hgenparticle_rapidity = new TH1F("gen_rapidity","Rapidity distribution", 41 , -20, 20 );

    Int_t nentries = (Int_t)mytree->GetEntries();

    for( Int_t i = 0 ; i < nentries ; i++)
    {
       mytree->GetEntry(i);
       hgenparticle_track->Fill(genparticle_track);
      for(int cnt = 0; cnt < genparticle_track ; cnt ++)
      {
       hgenparticle_charge->Fill(genparticle_charge[cnt]);
       hgenparticle_pdgId->Fill(genparticle_pdgId[cnt]);
       hgenparticle_px->Fill(genparticle_px[cnt]);
       hgenparticle_py->Fill(genparticle_py[cnt]);
       hgenparticle_pz->Fill(genparticle_pz[cnt]);
       hgenparticle_pt->Fill(genparticle_pt[cnt]);
       hgenparticle_mass->Fill(genparticle_mass[cnt]);
       hgenparticle_mt->Fill(genparticle_mt[cnt]);
       hgenparticle_phi->Fill(genparticle_phi[cnt]);
       hgenparticle_theta->Fill(genparticle_theta[cnt]);
       hgenparticle_eta->Fill(genparticle_eta[cnt]);
       hgenparticle_rapidity->Fill(genparticle_rapidity[cnt]);
      }
      /*
*/
    }
    hgenparticle_track->Draw();
    c1->SaveAs("ntuple.pdf(");
    hgenparticle_charge->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_pdgId->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_px->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_py->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_pz->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_pt->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_mass->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_mt->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_phi->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_theta->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_eta->Draw();
    c1->SaveAs("ntuple.pdf");
    hgenparticle_rapidity->Draw();
    c1->SaveAs("ntuple.pdf)");
}
void read_AK4PFJet()
{
//Branch for ak4pfjets
    //TFile *f = new TFile("/Users/chenzhenggang/Desktop/local_analysis/Ntuple_1_event.root");
    TFile *f = new TFile("/Users/chenzhenggang/Desktop/local_analysis/Ntuple.root");
    TTree *mytree = (TTree*)f->Get("Tree");
    TCanvas *c1=new TCanvas("c1","c1",1600,1600);
    
    Int_t ak4pfjets_track;
    Int_t ak4pfjets_charge[kMax];
    Int_t ak4pfjets_pdgId[kMax];
    Double32_t ak4pfjets_px[kMax];
    Double32_t ak4pfjets_py[kMax];
    Double32_t ak4pfjets_pz[kMax];
    Double32_t ak4pfjets_pt[kMax];
    Double32_t ak4pfjets_mass[kMax];
    Double32_t ak4pfjets_mt[kMax];
    Double32_t ak4pfjets_phi[kMax];
    Double32_t ak4pfjets_theta[kMax];
    Double32_t ak4pfjets_eta[kMax];
    Double32_t ak4pfjets_rapidity[kMax];
    
    mytree->SetBranchAddress("ak4pfjets_track",&ak4pfjets_track);
    mytree->SetBranchAddress("ak4pfjets_charge",&ak4pfjets_charge[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_pdgId",&ak4pfjets_pdgId);
    mytree->SetBranchAddress("ak4pfjets_px",&ak4pfjets_px[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_py",&ak4pfjets_py[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_pz",&ak4pfjets_pz[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_pt",&ak4pfjets_pt[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_mass",&ak4pfjets_mass[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_mt",&ak4pfjets_mt[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_phi",&ak4pfjets_phi[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_theta",&ak4pfjets_theta[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_eta",&ak4pfjets_eta[ak4pfjets_track]);
    mytree->SetBranchAddress("ak4pfjets_rapidity",&ak4pfjets_rapidity[ak4pfjets_track]);
    
    TH1F *hak4pfjets_track = new TH1F("ak4pfjets_track"," track distribution ", 31, 0 , 30);
    TH1F *hak4pfjets_charge = new TH1F("ak4pfjets_charge","charge distribution", 19, -8 , 10);
    TH1F *hak4pfjets_pdgId = new TH1F("ak4pfjets_pdgId","PDG_ID distribution", 5, -2 , 2);
    TH1F *hak4pfjets_px = new TH1F("ak4pfjets_px","Px distribution", 81, -4000, 4000 );
    TH1F *hak4pfjets_py = new TH1F("ak4pfjets_py","Py distribution", 81, -4000, 4000 );
    TH1F *hak4pfjets_pz = new TH1F("ak4pfjets_pz","Pz distribution", 120, -7000, 4000 );
    TH1F *hak4pfjets_pt = new TH1F("ak4pfjets_pt","Pt distribution", 41, 0, 4000 );
    TH1F *hak4pfjets_mass = new TH1F("ak4pfjets_mass","Mass distribution", 46, 0., 450. );
    TH1F *hak4pfjets_mt = new TH1F("ak4pfjets_mt","Transverse mass distribution", 41, 0, 4000 );
    TH1F *hak4pfjets_phi = new TH1F("ak4pfjets_phi","Phi distribution", 31, -TMath::Pi(),TMath::Pi() );
    TH1F *hak4pfjets_theta = new TH1F("ak4pfjets_theta","Theta distribution", 300, 0,  TMath::Pi());
    TH1F *hak4pfjets_eta = new TH1F("ak4pfjets_eta","Eta distribution", 201, -6., 6. );
    TH1F *hak4pfjets_rapidity = new TH1F("ak4pfjets_rapidity","Rapidity distribution",101 , -6., 6. );

    Int_t nentries = (Int_t)mytree->GetEntries();

    for( Int_t i = 0 ; i < nentries ; i++)
    {
       mytree->GetEntry(i);
       hak4pfjets_track->Fill(ak4pfjets_track);
      for(int cnt = 0; cnt < ak4pfjets_track ; cnt ++)
      {
       hak4pfjets_charge->Fill(ak4pfjets_charge[cnt]);
       hak4pfjets_pdgId->Fill(ak4pfjets_pdgId[cnt]);
       hak4pfjets_px->Fill(ak4pfjets_px[cnt]);
       hak4pfjets_py->Fill(ak4pfjets_py[cnt]);
       hak4pfjets_pz->Fill(ak4pfjets_pz[cnt]);
       hak4pfjets_pt->Fill(ak4pfjets_pt[cnt]);
       hak4pfjets_mass->Fill(ak4pfjets_mass[cnt]);
       hak4pfjets_mt->Fill(ak4pfjets_mt[cnt]);
       hak4pfjets_phi->Fill(ak4pfjets_phi[cnt]);
       hak4pfjets_theta->Fill(ak4pfjets_theta[cnt]);
       hak4pfjets_eta->Fill(ak4pfjets_eta[cnt]);
       hak4pfjets_rapidity->Fill(ak4pfjets_rapidity[cnt]);
      }
      /*
*/
    }
    hak4pfjets_track->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf(");
    hak4pfjets_charge->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_pdgId->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_px->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_py->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_pz->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_pt->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_mass->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_mt->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_phi->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_theta->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_eta->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf");
    hak4pfjets_rapidity->Draw();
    c1->SaveAs("ntuple_ak4pfjets.pdf)");
}


