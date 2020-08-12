#include <iostream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
void tree_testw();
void tree_testr();
void tree_testr2();
void tree_test();
void tree_test()
{
    tree_testw();
    tree_testr();
    tree_testr2();
}
void tree_testw()
{
    const int kMaxTrack = 1500;
    Int_t ntrack ;
    Int_t stat[kMaxTrack];
    Int_t sign[kMaxTrack];
    Double_t px[kMaxTrack];
    Double_t py[kMaxTrack];
    Double_t pz[kMaxTrack];
    Double_t pt[kMaxTrack];
    Double_t zv[kMaxTrack];
    Double_t chi2[kMaxTrack];
    Double_t sumstat;

    TFile f("tree_test.root","recreate");
    TTree *t = new TTree("t", "Reconst ntuple");
    t->Branch("ntrack", &ntrack, "ntrack/I");
    t->Branch("stat", stat, "stat[ntrack]/I");
    t->Branch("sign", sign, "sign[ntrack]/I");
    t->Branch("px", px, "px[ntrack]/D");
    t->Branch("py", py, "py[ntrack]/D");
    t->Branch("pz", pz, "pz[ntrack]/D");
    t->Branch("zv", zv, "zv[ntrack]/D");
    t->Branch("chi2",chi2,"chi2[ntrack]/D");
    
    TFile fr("tree_testf.root" , "recreate");
    TTree *tf = new TTree("tf" , "a friend tree");
    tf->Branch("ntrack" , &ntrack , "ntrack/I");
    tf->Branch("sumstat", &sumstat , " sumstat/D");
    tf->Branch("pt",pt,"pt[ntrack]/F");
    
    for(Int_t i = 0 ; i < 1000 ; i++)
    {
        Int_t nt = gRandom->Rndm() * (kMaxTrack - 1 );
        ntrack = nt ;
        sumstat = 0 ;
        for(Int_t n = 0 ; n < nt ; n++)
        {
            stat[ n ] = n%3;
            sign[ n ] = i%2;
            px[ n ] = gRandom->Gaus( 0 , 1 ) ;
            py[ n ] = gRandom->Gaus( 0 , 2 ) ;
            pz[ n ] = gRandom->Gaus(10 , 5 );
            zv[ n ] = gRandom->Gaus(100,2);
            chi2[ n ] = gRandom->Gaus( 0, .01);
            sumstat += chi2[ n ];
            pt[ n ] = TMath::Sqrt(px[ n ]* px[ n ] + py[ n ]*py[ n ]);
        }
        t->Fill();
        tf->Fill();

    }
    t->Print();
    f.cd();
    t->Write();
    fr.cd();
    tf->Write();

}

void tree_testr()
{
    TFile *f = new TFile("tree_test.root");
    TTree *t = (TTree * )f->Get("t");
    t->AddFriend("tf","tree_testf.root");
    t->Draw("pz","pt>3");
}
void tree_testr2()
{
    TPad *p =  new TPad("p", "p", 0.6 , 0.4 , 0.98 , 0.8);
    p->Draw();p->cd();
    TFile *f1 = new TFile("tree_test.root");
    TFile *f2 = new TFile("tree_testf.root");
    TTree *t = (TTree * )f1->Get("t");
    t->AddFriend("tf",f2);
    t->Draw("pz","pt>3");
}
