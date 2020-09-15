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
using namespace TMath;
const Int_t kMax = 10000;
void Find_min(Double_t [] , Int_t , Int_t , Double_t &);

void Delta_R_Distribution()
{
    TFile *f = new TFile("/Users/chenzhenggang/Desktop/local_analysis/Ntuple.root");
    TTree *mytree = (TTree*)f->Get("Tree");
    TCanvas *c1=new TCanvas("c1","c1",1600,1600);

    //Branch for genparticle
    Int_t genparticle_track;
    Int_t genparticle_charge[ kMax ];
    Int_t genparticle_pdgId[ kMax ];
    Int_t genparticle_status[ kMax ];
    Double_t genparticle_px[ kMax ];
    Double_t genparticle_py[ kMax ];
    Double_t genparticle_pz[ kMax ];
    Double_t genparticle_pt[ kMax ];
    Double_t genparticle_mass[ kMax ];
    Double_t genparticle_mt[ kMax ];
    Double_t genparticle_phi[ kMax ];
    Double_t genparticle_theta[ kMax ];
    Double_t genparticle_eta[ kMax ];
    Double_t genparticle_rapidity[ kMax ];

    mytree->SetBranchAddress("genparticle_track",&genparticle_track);
    mytree->SetBranchAddress("genparticle_charge",&genparticle_charge[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_pdgId",&genparticle_pdgId[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_status",&genparticle_status[ genparticle_track ]);
//    mytree->SetBranchAddress("genparticle_status",&genparticle_status[ genparticle_track ]) ;
    mytree->SetBranchAddress("genparticle_px",&genparticle_px[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_py",&genparticle_py[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_pz",&genparticle_pz[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_pt",&genparticle_pt[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_mass",&genparticle_mass[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_mt",&genparticle_mt[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_phi",&genparticle_phi[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_theta",&genparticle_theta[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_eta",&genparticle_eta[ genparticle_track ]);
    mytree->SetBranchAddress("genparticle_rapidity",&genparticle_rapidity[ genparticle_track ]);

    //Branch for AK4PFJets  
    
    Int_t ak4pfjets_track;
    Int_t ak4pfjets_charge[ kMax ];
    Int_t ak4pfjets_pdgId[ kMax ];
    Double_t ak4pfjets_px[ kMax ];
    Double_t ak4pfjets_py[ kMax ];
    Double_t ak4pfjets_pz[ kMax ];
    Double_t ak4pfjets_pt[ kMax ];
    Double_t ak4pfjets_mass[ kMax ];
    Double_t ak4pfjets_mt[ kMax ];
    Double_t ak4pfjets_phi[ kMax ];
    Double_t ak4pfjets_theta[ kMax ];
    Double_t ak4pfjets_eta[ kMax ];
    Double_t ak4pfjets_rapidity[ kMax ];
    
    mytree->SetBranchAddress("ak4pfjets_track", &ak4pfjets_track );
    mytree->SetBranchAddress("ak4pfjets_charge", &ak4pfjets_charge[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_pdgId", &ak4pfjets_pdgId );
    mytree->SetBranchAddress("ak4pfjets_px", &ak4pfjets_px[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_py", &ak4pfjets_py[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_pz", &ak4pfjets_pz[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_pt", &ak4pfjets_pt[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_mass", &ak4pfjets_mass[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_mt", &ak4pfjets_mt[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_phi", &ak4pfjets_phi[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_theta", &ak4pfjets_theta[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_eta", &ak4pfjets_eta[ ak4pfjets_track ] );
    mytree->SetBranchAddress("ak4pfjets_rapidity", &ak4pfjets_rapidity[ ak4pfjets_track ] );
    
    TH1D *Minimum_Delta_R_Distribution = new TH1D("Mini_Delta_R","Minimum Delta R Distribution ", 100, 0 , 1);
    TH1D *Delta_R_Distribution = new TH1D("Delta_R","Delta R Distribution ", 100, 0 , 1);
    Int_t nentries = (Int_t)mytree->GetEntries();
    for( Int_t ev = 0 ; ev < nentries ; ev++)
    {
       mytree->GetEntry(ev);
       int ctr_gen = 0;
      
      for(int i = 0; i < genparticle_track ; i ++)
      {
      
       if(genparticle_pdgId[ i ] < 4 and genparticle_pdgId[ i ] > 0 and genparticle_status[ i ] == 23 )
       {
           Double_t dR[ ak4pfjets_track ] ;
           Double_t minimum_DeltaR=100000000 ; 
           for(Int_t j = 0; j < ak4pfjets_track ; j ++)
                  {
                        dR[ j ] = Sqrt(Power((ak4pfjets_eta[ j ] - genparticle_eta[ i ]),2) + Power((ak4pfjets_phi[ j ] - genparticle_phi[ i ]), 2 ) );
                         Minimum_Delta_R_Distribution->Fill( dR[ j ] );
                  }

           Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR );
           Minimum_Delta_R_Distribution->Fill( minimum_DeltaR );
       }
      }
    }
    c1->cd(); 
    Minimum_Delta_R_Distribution->Draw();
    c1->SaveAs("Delta_R(.pdf");
    Delta_R_Distribution->Draw();
    c1->SaveAs("Delta_R).pdf");
    
    c1->Close();
    }

void Find_min(Double_t arr[] , const Int_t size , Int_t pos ,Double_t &min )
{
    if(pos != size - 1)
    {
        Find_min( arr, size , pos + 1 , min );
    }
    if( min > arr[ pos ])
        min = arr[ pos ];
}

