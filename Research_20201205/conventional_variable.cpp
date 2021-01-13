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
const double hgcal_upper_eta = 3.0;

const double part_dR_cut =0.1;
const Double_t R_cut = 0.4;
const Double_t constraint_coef_pt = .8;
const double clornot = 0.03;

const double canvas_width = 900;
const double font_size = 0.055;
const double font_offset= 0.9;
const double legend_size=0.025;

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
double SC_Energy=0.0;//selection cut for clsuters energy 
double SC_neutral_part_Pt = 1.;

void conventional_variable()
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
    std::vector<Double_t> *PFConstituents_pt;
    std::vector<Int_t> *PFConstituents_charge;
    std::vector<Int_t> *PFConstituents_pdgID;
    mytree->SetBranchAddress("PFConstituents_eta", &PFConstituents_eta,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_phi", &PFConstituents_phi,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_charge", &PFConstituents_charge,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_pt", &PFConstituents_pt,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_E", &PFConstituents_E,&PFConstituents);
    mytree->SetBranchAddress("PFConstituents_pdgId",&PFConstituents_pdgID,&PFConstituents);
    
    Int_t nentries = (Int_t)mytree->GetEntries();
    //Setting for historgram 
    TH1D *Quark_charged_multi;
    TH1D *Gluon_charged_multi;
    TH1D *Quark_charged_multi_CL;
    TH1D *Gluon_charged_multi_CL;
    
    TH1D *Quark_frag;
    TH1D *Gluon_frag;
    TH1D *Quark_frag_CL;
    TH1D *Gluon_frag_CL;
    
    TH1D *Quark_major_axis;
    TH1D *Gluon_major_axis;
    TH1D *Quark_major_axis_CL;
    TH1D *Gluon_major_axis_CL;
    
    TH1D *Quark_minor_axis;
    TH1D *Gluon_minor_axis;
    TH1D *Quark_minor_axis_CL;
    TH1D *Gluon_minor_axis_CL;
    
    TLegend legend_multi_CL;
    TLegend legend_multi;
    TLegend legend_frag_CL;
    TLegend legend_frag;
    TLegend legend_major_axis_CL;
    TLegend legend_major_axis;
    TLegend legend_minor_axis_CL;
    TLegend legend_minor_axis;

    char titleg_multi_CL[80];
    char titleq_multi_CL[80];
    char nameq_multi_CL[80];
    char nameg_multi_CL[80];
    
    char titleg_multi[80];
    char titleq_multi[80];
    char nameq_multi[80];
    char nameg_multi[80];
    
    sprintf(titleg_multi,"Charged Particles Number");
    sprintf(titleq_multi,"Charged Particles Number");
    
    sprintf(nameg_multi,"Charged Particles Number for quark");
    sprintf(nameq_multi,"Charged Particles Number for gluon");
  
    sprintf(titleg_multi_CL,"Charged Clusters Number");
    sprintf(titleq_multi_CL,"Charged Clusters Number");
    
    sprintf(nameg_multi_CL,"Charged Clusters Number for quark");
    sprintf(nameq_multi_CL,"Charged Clusters Number for gluon");

    char titleg_frag_CL[80];
    char titleq_frag_CL[80];
    char nameq_frag_CL[80];
    char nameg_frag_CL[80];
    
    char titleg_frag[80];
    char titleq_frag[80];
    char nameq_frag[80];
    char nameg_frag[80];
    
    sprintf(titleg_frag,"Particles:Fragmentation Function");
    sprintf(titleq_frag,"Particles:Fragmentation Function");
    
    sprintf(nameg_frag,"Particles:Fragmentation Function for quark");
    sprintf(nameq_frag,"Particles:Fragmentation Function for gluon");
  
    sprintf(titleg_frag_CL,"Clusters:Fragmentation Function");
    sprintf(titleq_frag_CL,"Clusters:Fragmentation Function");
    
    sprintf(nameg_frag_CL,"Clusters:Fragmentation Function for quark");
    sprintf(nameq_frag_CL,"Clusters:Fragmentation Function for gluon");
    
    char titleg_major_axis_CL[80];
    char titleq_major_axis_CL[80];
    char nameq_major_axis_CL[80];
    char nameg_major_axis_CL[80];
    
    char titleg_major_axis[80];
    char titleq_major_axis[80];
    char nameq_major_axis[80];
    char nameg_major_axis[80];
    
    sprintf(titleg_major_axis,"Particles: Jet Major axis");
    sprintf(titleq_major_axis,"Particles: Jet Major axis");
    
    sprintf(nameg_major_axis,"Particles: Jet Major axis for quark");
    sprintf(nameq_major_axis,"Particles: Jet Major axis for gluon");
  
    sprintf(titleg_major_axis_CL,"Clusters: Jet Major axis");
    sprintf(titleq_major_axis_CL,"Clusters: Jet Major axis");
    
    sprintf(nameg_major_axis_CL,"Clusters: Jet Major axis for quark");
    sprintf(nameq_major_axis_CL,"Clusters: Jet Major axis for gluon");
    
    char titleg_minor_axis_CL[80];
    char titleq_minor_axis_CL[80];
    char nameq_minor_axis_CL[80];
    char nameg_minor_axis_CL[80];
    
    char titleg_minor_axis[80];
    char titleq_minor_axis[80];
    char nameq_minor_axis[80];
    char nameg_minor_axis[80];
    
    sprintf(titleg_minor_axis,"Particles: Jet Minor axis");
    sprintf(titleq_minor_axis,"Particles: Jet Minor axis");
    
    sprintf(nameg_minor_axis,"Particles: Jet Minor axis for quark");
    sprintf(nameq_minor_axis,"Particles: Jet Minor axis for gluon");
  
    sprintf(titleg_minor_axis_CL,"Clusters: Jet Minor axis");
    sprintf(titleq_minor_axis_CL,"Clusters: Jet Minor axis");
    
    sprintf(nameg_minor_axis_CL,"Clusters: Jet Minor axis for quark");
    sprintf(nameq_minor_axis_CL,"Clusters: Jet Minor axis for gluon");
    
    Quark_charge_multi=new TH1D(nameq_multi,titleq_multi,30,0.,60.);
    Quark_charge_multi->SetFillColorAlpha(kRed,0.15);
    Quark_charge_multi->GetYaxis()->SetTitle("Probability");
    Quark_charge_multi->GetXaxis()->SetTitle("Numbers");
    Quark_charge_multi->GetXaxis()->SetTitleSize(font_size);
    Quark_charge_multi->GetYaxis()->SetTitleSize(font_size);
    Quark_charge_multi->GetXaxis()->SetTitleOffset(font_offset);
    Quark_charge_multi->GetYaxis()->SetTitleOffset(font_offset);
    Quark_charge_multi->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_charge_multi->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_charge_multi->SetLineWidth(2);
    Quark_charge_multi->SetLineColor(kRed);

    Gluon_charge_multi=new TH1D(nameg_multi,titleg_multi,30,0.,60.);
    Gluon_charge_multi->SetFillColorAlpha(kBlue,0.15);
    Gluon_charge_multi->GetYaxis()->SetTitle("Probability");
    Gluon_charge_multi->GetXaxis()->SetTitle("Numbers");
    Gluon_charge_multi->GetXaxis()->SetTitleSize(font_size);
    Gluon_charge_multi->GetYaxis()->SetTitleSize(font_size);
    Gluon_charge_multi->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_charge_multi->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_charge_multi->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_charge_multi->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_charge_multi->SetLineWidth(2);
    Gluon_charge_multi->SetLineColor(kBlue);

   
    legend_multi.SetX1NDC(0.35);
    legend_multi.SetY1NDC(0.7);
    legend_multi.SetX2NDC(0.88);
    legend_multi.SetY2NDC(0.88);
    legend_multi.SetHeader("Representation");
    legend_multi.SetTextSize(legend_size);
    legend_multi.AddEntry(Quark_charge_multi,nameg_multi,"f");
    legend_multi.AddEntry(Gluon_charge_multi,nameq_multi,"f");
    
    Quark_charge_multi_CL=new TH1D(nameq_multi_CL,titleq_multi_CL,15,0.,15000.);
    Quark_charge_multi_CL->SetFillColorAlpha(kRed,0.15);
    Quark_charge_multi_CL->GetYaxis()->SetTitle("Probability");
    Quark_charge_multi_CL->GetXaxis()->SetTitle("Numbers");
    Quark_charge_multi_CL->GetXaxis()->SetTitleSize(font_size);
    Quark_charge_multi_CL->GetYaxis()->SetTitleSize(font_size);
    Quark_charge_multi_CL->GetXaxis()->SetTitleOffset(font_offset);
    Quark_charge_multi_CL->GetYaxis()->SetTitleOffset(font_offset);
    Quark_charge_multi_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_charge_multi_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_charge_multi_CL->SetLineWidth(2);
    Quark_charge_multi_CL->SetLineColor(kRed);

    Gluon_charge_multi_CL=new TH1D(nameg_multi_CL,titleg_multi_CL,15,0.,15000.);
    Gluon_charge_multi_CL->SetFillColorAlpha(kBlue,0.15);
    Gluon_charge_multi_CL->GetYaxis()->SetTitle("Probability");
    Gluon_charge_multi_CL->GetXaxis()->SetTitle("Numbers");
    Gluon_charge_multi_CL->GetXaxis()->SetTitleSize(font_size);
    Gluon_charge_multi_CL->GetYaxis()->SetTitleSize(font_size);
    Gluon_charge_multi_CL->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_charge_multi_CL->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_charge_multi_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_charge_multi_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_charge_multi_CL->SetLineWidth(2);
    Gluon_charge_multi_CL->SetLineColor(kBlue);

   
    legend_multi_CL.SetX1NDC(0.35);
    legend_multi_CL.SetY1NDC(0.7);
    legend_multi_CL.SetX2NDC(0.88);
    legend_multi_CL.SetY2NDC(0.88);
    legend_multi_CL.SetHeader("Representation");
    legend_multi_CL.SetTextSize(legend_size);
    legend_multi_CL.AddEntry(Quark_charge_multi_CL,nameg_multi_CL,"f");
    legend_multi_CL.AddEntry(Gluon_charge_multi_CL,nameq_multi_CL,"f");

    Quark_frag=new TH1D(nameq_multi,titleq_multi,50,0.,1.);
    Quark_frag->SetFillColorAlpha(kRed,0.15);
    Quark_frag->GetYaxis()->SetTitle("Probability");
    Quark_frag->GetXaxis()->SetTitle("P_{f}");
    Quark_frag->GetXaxis()->SetTitleSize(font_size);
    Quark_frag->GetYaxis()->SetTitleSize(font_size);
    Quark_frag->GetXaxis()->SetTitleOffset(font_offset);
    Quark_frag->GetYaxis()->SetTitleOffset(font_offset);
    Quark_frag->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_frag->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_frag->SetLineWidth(2);
    Quark_frag->SetLineColor(kRed);

    Gluon_frag=new TH1D(nameg_multi,titleg_multi,50,0.,1.);
    Gluon_frag->SetFillColorAlpha(kBlue,0.15);
    Gluon_frag->GetYaxis()->SetTitle("Probability");
    Gluon_frag->GetXaxis()->SetTitle("P_{f}");
    Gluon_frag->GetXaxis()->SetTitleSize(font_size);
    Gluon_frag->GetYaxis()->SetTitleSize(font_size);
    Gluon_frag->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_frag->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_frag->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_frag->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_frag->SetLineWidth(2);
    Gluon_frag->SetLineColor(kBlue);

   
    legend_frag.SetX1NDC(0.35);
    legend_frag.SetY1NDC(0.7);
    legend_frag.SetX2NDC(0.88);
    legend_frag.SetY2NDC(0.88);
    legend_frag.SetHeader("Representation");
    legend_frag.SetTextSize(legend_size);
    legend_frag.AddEntry(Quark_charge_multi,nameg_multi,"f");
    legend_frag.AddEntry(Gluon_charge_multi,nameq_multi,"f");
    
    Quark_frag_CL=new TH1D(nameq_multi_CL,titleq_multi_CL,50,0.,0.1);
    Quark_frag_CL->SetFillColorAlpha(kRed,0.15);
    Quark_frag_CL->GetYaxis()->SetTitle("Probability");
    Quark_frag_CL->GetXaxis()->SetTitle("P_{f}");
    Quark_frag_CL->GetXaxis()->SetTitleSize(font_size);
    Quark_frag_CL->GetYaxis()->SetTitleSize(font_size);
    Quark_frag_CL->GetXaxis()->SetTitleOffset(font_offset);
    Quark_frag_CL->GetYaxis()->SetTitleOffset(font_offset);
    Quark_frag_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_frag_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_frag_CL->SetLineWidth(2);
    Quark_frag_CL->SetLineColor(kRed);

    Gluon_frag_CL=new TH1D(nameg_multi_CL,titleg_multi_CL,50,0.,0.1);
    Gluon_frag_CL->SetFillColorAlpha(kBlue,0.15);
    Gluon_frag_CL->GetYaxis()->SetTitle("Probability");
    Gluon_frag_CL->GetXaxis()->SetTitle("P_{f}");
    Gluon_frag_CL->GetXaxis()->SetTitleSize(font_size);
    Gluon_frag_CL->GetYaxis()->SetTitleSize(font_size);
    Gluon_frag_CL->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_frag_CL->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_frag_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_frag_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_frag_CL->SetLineWidth(2);
    Gluon_frag_CL->SetLineColor(kBlue);

    legend_frag_CL.SetX1NDC(0.35);
    legend_frag_CL.SetY1NDC(0.7);
    legend_frag_CL.SetX2NDC(0.88);
    legend_frag_CL.SetY2NDC(0.88);
    legend_frag_CL.SetHeader("Representation");
    legend_frag_CL.SetTextSize(legend_size);
    legend_frag_CL.AddEntry(Quark_charge_multi_CL,nameg_multi_CL,"f");
    legend_frag_CL.AddEntry(Gluon_charge_multi_CL,nameq_multi_CL,"f");

    Quark_major_axis=new TH1D(nameq_major_axis,titleq_major_axis,30,0.,0.3);
    Quark_major_axis->SetFillColorAlpha(kRed,0.15);
    Quark_major_axis->GetYaxis()->SetTitle("Probabilitiy");
    Quark_major_axis->GetXaxis()->SetTitle("");
    Quark_major_axis->GetXaxis()->SetTitleSize(font_size);
    Quark_major_axis->GetYaxis()->SetTitleSize(font_size);
    Quark_major_axis->GetXaxis()->SetTitleOffset(font_offset);
    Quark_major_axis->GetYaxis()->SetTitleOffset(font_offset);
    Quark_major_axis->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_major_axis->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_major_axis->SetLineWidth(2);
    Quark_major_axis->SetLineColor(kRed);

    Gluon_major_axis=new TH1D(nameg_major_axis,titleg_major_axis,30,0.,0.3);
    Gluon_major_axis->SetFillColorAlpha(kBlue,0.15);
    Gluon_major_axis->GetYaxis()->SetTitle("Probability");
    Gluon_major_axis->GetXaxis()->SetTitle("");
    Gluon_major_axis->GetXaxis()->SetTitleSize(font_size);
    Gluon_major_axis->GetYaxis()->SetTitleSize(font_size);
    Gluon_major_axis->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_major_axis->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_major_axis->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_major_axis->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_major_axis->SetLineWidth(2);
    Gluon_major_axis->SetLineColor(kBlue);

   
    legend_major_axis.SetX1NDC(0.35);
    legend_major_axis.SetY1NDC(0.7);
    legend_major_axis.SetX2NDC(0.88);
    legend_major_axis.SetY2NDC(0.88);
    legend_major_axis.SetHeader("Representation");
    legend_major_axis.SetTextSize(legend_size);
    legend_major_axis.AddEntry(Quark_major_axis,nameg_major_axis,"f");
    legend_major_axis.AddEntry(Gluon_major_axis,nameq_major_axis,"f");
    
    Quark_major_axis_CL=new TH1D(nameq_major_axis_CL,titleq_major_axis_CL,50,0.,0.05);
    Quark_major_axis_CL->SetFillColorAlpha(kRed,0.15);
    Quark_major_axis_CL->GetYaxis()->SetTitle("Probability");
    Quark_major_axis_CL->GetXaxis()->SetTitle("");
    Quark_major_axis_CL->GetXaxis()->SetTitleSize(font_size);
    Quark_major_axis_CL->GetYaxis()->SetTitleSize(font_size);
    Quark_major_axis_CL->GetXaxis()->SetTitleOffset(font_offset);
    Quark_major_axis_CL->GetYaxis()->SetTitleOffset(font_offset);
    Quark_major_axis_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_major_axis_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_major_axis_CL->SetLineWidth(2);
    Quark_major_axis_CL->SetLineColor(kRed);

    Gluon_major_axis_CL=new TH1D(nameg_major_axis_CL,titleg_major_axis_CL,50,0.,0.05);
    Gluon_major_axis_CL->SetFillColorAlpha(kBlue,0.15);
    Gluon_major_axis_CL->GetYaxis()->SetTitle("Probability");
    Gluon_major_axis_CL->GetXaxis()->SetTitle("");
    Gluon_major_axis_CL->GetXaxis()->SetTitleSize(font_size);
    Gluon_major_axis_CL->GetYaxis()->SetTitleSize(font_size);
    Gluon_major_axis_CL->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_major_axis_CL->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_major_axis_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_major_axis_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_major_axis_CL->SetLineWidth(2);
    Gluon_major_axis_CL->SetLineColor(kBlue);

   
    legend_major_axis_CL.SetX1NDC(0.35);
    legend_major_axis_CL.SetY1NDC(0.7);
    legend_major_axis_CL.SetX2NDC(0.88);
    legend_major_axis_CL.SetY2NDC(0.88);
    legend_major_axis_CL.SetHeader("Representation");
    legend_major_axis_CL.SetTextSize(legend_size);
    legend_major_axis_CL.AddEntry(Quark_major_axis_CL,nameg_major_axis_CL,"f");
    legend_major_axis_CL.AddEntry(Gluon_major_axis_CL,nameq_major_axis_CL,"f");
    
    Quark_minor_axis=new TH1D(nameq_minor_axis,titleq_minor_axis,50,0.,0.1);
    Quark_minor_axis->SetFillColorAlpha(kRed,0.15);
    Quark_minor_axis->GetYaxis()->SetTitle("Probabilitiy");
    Quark_minor_axis->GetXaxis()->SetTitle("");
    Quark_minor_axis->GetXaxis()->SetTitleSize(font_size);
    Quark_minor_axis->GetYaxis()->SetTitleSize(font_size);
    Quark_minor_axis->GetXaxis()->SetTitleOffset(font_offset);
    Quark_minor_axis->GetYaxis()->SetTitleOffset(font_offset);
    Quark_minor_axis->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_minor_axis->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_minor_axis->SetLineWidth(2);
    Quark_minor_axis->SetLineColor(kRed);

    Gluon_minor_axis=new TH1D(nameg_minor_axis,titleg_minor_axis,50,0.,0.1);
    Gluon_minor_axis->SetFillColorAlpha(kBlue,0.15);
    Gluon_minor_axis->GetYaxis()->SetTitle("Probability");
    Gluon_minor_axis->GetXaxis()->SetTitle("");
    Gluon_minor_axis->GetXaxis()->SetTitleSize(font_size);
    Gluon_minor_axis->GetYaxis()->SetTitleSize(font_size);
    Gluon_minor_axis->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_minor_axis->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_minor_axis->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_minor_axis->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_minor_axis->SetLineWidth(2);
    Gluon_minor_axis->SetLineColor(kBlue);

   
    legend_minor_axis.SetX1NDC(0.35);
    legend_minor_axis.SetY1NDC(0.7);
    legend_minor_axis.SetX2NDC(0.88);
    legend_minor_axis.SetY2NDC(0.88);
    legend_minor_axis.SetHeader("Representation");
    legend_minor_axis.SetTextSize(legend_size);
    legend_minor_axis.AddEntry(Quark_minor_axis,nameg_minor_axis,"f");
    legend_minor_axis.AddEntry(Gluon_minor_axis,nameq_minor_axis,"f");
    
    Quark_minor_axis_CL=new TH1D(nameq_minor_axis_CL,titleq_minor_axis_CL,50,0.,0.005);
    Quark_minor_axis_CL->SetFillColorAlpha(kRed,0.15);
    Quark_minor_axis_CL->GetYaxis()->SetTitle("Probability");
    Quark_minor_axis_CL->GetXaxis()->SetTitle("");
    Quark_minor_axis_CL->GetXaxis()->SetTitleSize(font_size);
    Quark_minor_axis_CL->GetYaxis()->SetTitleSize(font_size);
    Quark_minor_axis_CL->GetXaxis()->SetTitleOffset(font_offset);
    Quark_minor_axis_CL->GetYaxis()->SetTitleOffset(font_offset);
    Quark_minor_axis_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Quark_minor_axis_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Quark_minor_axis_CL->SetLineWidth(2);
    Quark_minor_axis_CL->SetLineColor(kRed);

    Gluon_minor_axis_CL=new TH1D(nameg_minor_axis_CL,titleg_minor_axis_CL,50,0.,0.005);
    Gluon_minor_axis_CL->SetFillColorAlpha(kBlue,0.15);
    Gluon_minor_axis_CL->GetYaxis()->SetTitle("Probability");
    Gluon_minor_axis_CL->GetXaxis()->SetTitle("");
    Gluon_minor_axis_CL->GetXaxis()->SetTitleSize(font_size);
    Gluon_minor_axis_CL->GetYaxis()->SetTitleSize(font_size);
    Gluon_minor_axis_CL->GetXaxis()->SetTitleOffset(font_offset);
    Gluon_minor_axis_CL->GetYaxis()->SetTitleOffset(font_offset);
    Gluon_minor_axis_CL->GetXaxis()->SetLabelSize(0.4*font_size);
    Gluon_minor_axis_CL->GetYaxis()->SetLabelSize(0.4*font_size);
    Gluon_minor_axis_CL->SetLineWidth(2);
    Gluon_minor_axis_CL->SetLineColor(kBlue);

   
    legend_minor_axis_CL.SetX1NDC(0.35);
    legend_minor_axis_CL.SetY1NDC(0.7);
    legend_minor_axis_CL.SetX2NDC(0.88);
    legend_minor_axis_CL.SetY2NDC(0.88);
    legend_minor_axis_CL.SetHeader("Representation");
    legend_minor_axis_CL.SetTextSize(legend_size);
    legend_minor_axis_CL.AddEntry(Quark_minor_axis_CL,nameg_minor_axis_CL,"f");
    legend_minor_axis_CL.AddEntry(Gluon_minor_axis_CL,nameq_minor_axis_CL,"f");
    
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
                       dR[j] = TMath::Sqrt(TMath::Power((ak4pfjets_eta[j] - genparticle_eta[i]),2) + \
                               TMath::Power((ak4pfjets_phi[j] - genparticle_phi[i]), 2 ) );
                   }
                   
                   Find_min( dR, ak4pfjets_track , 0 , minimum_DeltaR , min_pos , flag );
                   if( minimum_DeltaR < 0.1 && TMath::Abs(ak4pfjets_eta[min_pos]) < hgcal_upper_eta \
                           && TMath::Abs( ak4pfjets_eta[min_pos])>hgcal_lower_eta  \
                           && ak4pfjets_pt[min_pos]>genparticle_pt[i]*constraint_coef_pt\
                           )
                   {
                       //Genparticle is Gluon or Quark
                       //Random choose one quark jet or gluon jet from events
                       int charged_cluster = 0;
                       int charged_particle = 0;
                       double frag_particle = 0.;
                       double frag_cluster = 0.;
                       double E_total_particle = 0.;
                       double E_total_cluster = 0.;
                       double width_particle = 0. ;
                       double width_cluster = 0. ;
                       
                       double M11_particle = 0.;
                       double M22_particle = 0.;
                       double Moff_particle = 0.;
                       
                       double M11_cluster = 0.;
                       double M22_cluster = 0.;
                       double Moff_cluster = 0.;
                       
                       double width_minor_particle = 0.;
                       double width_major_particle = 0.;
                       double width_minor_cluster = 0.;
                       double width_major_cluster = 0.;
                       double average_width_particle = 0.;
                       double average_width_cluster = 0.;
                       int number_particle = 0; 
                       int number_cluster = 0;
                       double avg_eta_particle = 0.;
                       double avg_phi_particle = 0.;
                       double E2_particle = 0.; 
                       double avg_eta_cluster = 0.;
                       double avg_phi_cluster = 0.;
                       double E2_cluster = 0.; 
                       Int_t particle_index;
                       Int_t index_initial;
                       if(min_pos >0 )
                       {
                           particle_index = ak4pfjets_tagging->at(min_pos-1);
                       }
                       else
                       {
                           particle_index=0;
                       }
                       index_initial = particle_index;
                       Int_t particle_index_end = ak4pfjets_tagging->at(min_pos);
                   //Loop for particles 
                       for(;particle_index<particle_index_end;particle_index++)
                       {
                           double energy = PFConstituents_E->at(particle_index);
                           double eta = PFConstituents_eta->at(particle_index);
                           double phi = PFConstituents_phi->at(particle_index);
                           double pt = PFConstituents_pt->at(particle_index);
                           if(PFConstituents_charge->at(particle_index)!=0)
                           {
                               charged_particle++;
                               number_particle++;
                               E_total_particle += energy;
                               frag_particle += energy*energy;
                               avg_phi_particle += phi;
                               avg_eta_particle += eta;
                               E2_particle += energy*energy;
                           }
                           else if(pt > SC_neutral_part_Pt)
                           {
                               number_particle++;
                               E_total_particle += energy;
                               frag_particle += energy*energy;
                               avg_phi_particle += phi;
                               avg_eta_particle += eta;
                               E2_particle += energy*energy;
                           }
                           for(Int_t k = 0 ; k < LayerCluster_track; k++ )
                           {
                               Double_t dR_=500;
                               int min_index=-1;
                               double deltaR=0;
                               for(int layer=0;layer<hgcallayer_max;layer++)
                               {
                                   if(LayerCluster_layer[k]==layer+1)
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
                                           if(deltaR < clornot)
                                           {
                                               if((electron||charged_Had)&&LayerCluster_E[k]>SC_Energy)
                                               {
                                                   charged_cluster++;
                                                   E_total_cluster += LayerCluster_E[k];
                                                   frag_cluster += LayerCluster_E[k]*LayerCluster_E[k];
                                                   avg_phi_cluster += phi;
                                                   avg_eta_cluster += eta;
                                                   number_cluster++;
                                                   E2_cluster += LayerCluster_E[k]*LayerCluster_E[k];
                                                   LayerCluster_tagger.at(k)=true;
                                               }
                                               else if(LayerCluster_E[k] > SC_Energy)
                                               {
                                                   E_total_cluster += LayerCluster_E[k];
                                                   frag_cluster += LayerCluster_E[k]*LayerCluster_E[k];
                                                   avg_phi_cluster += phi;
                                                   avg_eta_cluster += eta;
                                                   number_cluster++;
                                                   E2_cluster = LayerCluster_E[k]*LayerCluster_E[k];
                                                   LayerCluster_tagger.at(k)=true;
                                               }
                                           }
                                       }

                                   }
                               }
                           }
                       }
                       avg_phi_particle /= number_particle;
                       avg_eta_particle /= number_particle;
                       avg_phi_cluster /= number_cluster;
                       avg_eta_cluster /= number_cluster;
                       particle_index = index_initial;
                       for(;particle_index<particle_index_end;particle_index++)
                       {
                           double energy = PFConstituents_E->at(particle_index);
                           double pt = PFConstituents_pt->at(particle_index);
                           double deta = PFConstituents_eta->at(particle_index)-avg_eta_particle;
                           double dphi = PFConstituents_phi->at(particle_index)-avg_phi_particle;
                           if(PFConstituents_charge->at(particle_index)!=0)
                           {
                               M11_particle += energy*energy*deta*deta;
                               M22_particle += energy*energy*dphi*dphi;
                               Moff_particle += -energy*energy*dphi*deta;
                           }
                           else
                           {
                               if(pt > SC_neutral_part_Pt)
                               {
                                   M11_particle += energy*energy*deta*deta;
                                   M22_particle += energy*energy*dphi*dphi;
                                   Moff_particle += -energy*energy*dphi*deta;
                               }
                           }
                       }
                       for(int  k = 0 ; k < LayerCluster_track ; k++)
                       {
                           if( LayerCluster_tagger.at(k))
                           {
                               double energy = LayerCluster_E[k];
                               double deta = LayerCluster_eta[k]-avg_eta_cluster;
                               double dphi = LayerCluster_phi[k]-avg_phi_cluster;
                               M11_cluster += energy*energy*deta*deta; 
                               M22_cluster += energy*energy*dphi*dphi;
                               Moff_cluster += -energy*energy*dphi*deta;
                           }
                       }
                       
                       frag_particle = TMath::Sqrt(frag_particle);
                       frag_particle /= E_total_particle;
                       
                       frag_cluster = TMath::Sqrt(frag_cluster);
                       frag_cluster /= E_total_cluster;
                       
                       width_minor_particle = ((M11_particle+M22_particle)-TMath::Sqrt((M11_particle+M22_particle)*(M11_particle+M22_particle)\
                                   -4*(M11_particle*M22_particle-Moff_particle*Moff_particle)))/2.;
                       width_major_particle = ((M11_particle+M22_particle)+TMath::Sqrt((M11_particle+M22_particle)*(M11_particle+M22_particle)\
                                   -4*(M11_particle*M22_particle-Moff_particle*Moff_particle)))/2.;
                       width_minor_particle /= E2_particle;
                       width_major_particle /= E2_particle;
                       width_minor_particle = TMath::Sqrt(width_minor_particle);
                       width_major_particle = TMath::Sqrt(width_major_particle);
                       
                       width_minor_cluster = ((M11_cluster+M22_cluster)-TMath::Sqrt((M11_cluster+M22_cluster)*(M11_cluster+M22_cluster)\
                                   -4*(M11_cluster*M22_cluster-Moff_cluster*Moff_cluster)))/2.;
                       width_major_cluster = ((M11_cluster+M22_cluster)+TMath::Sqrt((M11_cluster+M22_cluster)*(M11_cluster+M22_cluster)\
                                   -4*(M11_cluster*M22_cluster-Moff_cluster*Moff_cluster)))/2.;
                       width_minor_cluster /= E2_cluster;
                       width_major_cluster /= E2_cluster;
                       width_minor_cluster = TMath::Sqrt(width_minor_cluster);
                       width_major_cluster = TMath::Sqrt(width_major_cluster);
                       if(quark)
                       {
                           Quark_charge_multi_CL->Fill(charged_cluster);
                           Quark_charge_multi->Fill(charged_particle);
                           Quark_major_axis_CL->Fill(width_major_cluster);
                           Quark_major_axis->Fill(width_major_particle);
                           Quark_minor_axis_CL->Fill(width_minor_cluster);
                           Quark_minor_axis->Fill(width_minor_particle);
                           Quark_frag_CL->Fill(frag_cluster);
                           Quark_frag->Fill(frag_particle);
                       }
                       else if(gluon)
                       {
                           Gluon_charge_multi_CL->Fill(charged_cluster);
                           Gluon_charge_multi->Fill(charged_particle);
                           Gluon_major_axis_CL->Fill(width_major_cluster);
                           Gluon_major_axis->Fill(width_major_particle);
                           Gluon_minor_axis_CL->Fill(width_minor_cluster);
                           Gluon_minor_axis->Fill(width_minor_particle);
                           Gluon_frag_CL->Fill(frag_cluster);
                           Gluon_frag->Fill(frag_particle);
                       }
                  }
               }
          }
   }
    //Multiplicity
    c1->cd();
    char q_title_multi_CL[80] ="";
    sprintf(q_title_multi_CL,"Multiplicity.pdf(");
    double scale_q_CL_multi = 1./Quark_charge_multi_CL->Integral();
    double scale_g_CL_multi = 1./Gluon_charge_multi_CL->Integral();
        
    Quark_charge_multi_CL->Scale(scale_q_CL_multi);
    Gluon_charge_multi_CL->Scale(scale_g_CL_multi);
    Gluon_charge_multi_CL->SetMaximum(2.5*Quark_charge_multi_CL->GetMaximum());
    Gluon_charge_multi_CL->Draw();
    Quark_charge_multi_CL->Draw("SAME");
    legend_multi_CL.Draw("SAME");
        
    c1->SaveAs(q_title_multi_CL);
    c1->Clear();
    
    char q_title_multi[80] ="";
    sprintf(q_title_multi,"Multiplicity.pdf)");
    double scale_q_part_multi = 1./Quark_charge_multi->Integral();
    double scale_g_part_multi = 1./Gluon_charge_multi->Integral();
        
    Quark_charge_multi->Scale(scale_q_part_multi);
    Gluon_charge_multi->Scale(scale_g_part_multi);
    Gluon_charge_multi->SetMaximum(2.5*Quark_charge_multi->GetMaximum());
    Gluon_charge_multi->Draw();
    Quark_charge_multi->Draw("SAME");
    legend_multi.Draw("SAME");
        
    c1->SaveAs(q_title_multi);
    c1->Clear();
    //fragmentation
    char q_title_frag_CL[80] ="";
    sprintf(q_title_frag_CL,"Fragmentation.pdf(");
    double scale_q_CL_frag = 1./Quark_frag_CL->Integral();
    double scale_g_CL_frag = 1./Gluon_frag_CL->Integral();
        
    Quark_frag_CL->Scale(scale_q_CL_frag);
    Gluon_frag_CL->Scale(scale_g_CL_frag);
    Gluon_frag_CL->SetMaximum(2.5*Quark_frag_CL->GetMaximum());
    Gluon_frag_CL->Draw();
    Quark_frag_CL->Draw("SAME");
    legend_frag_CL.Draw("SAME");
        
    c1->SaveAs(q_title_frag_CL);
    c1->Clear();
    
    char q_title_frag[80] ="";
    sprintf(q_title_frag,"Fragmentation.pdf)");
    double scale_q_part_frag = 1./Quark_frag->Integral();
    double scale_g_part_frag = 1./Gluon_frag->Integral();
        
    Quark_frag->Scale(scale_q_part_frag);
    Gluon_frag->Scale(scale_g_part_frag);
    Gluon_frag->SetMaximum(2.5*Quark_frag->GetMaximum());
    Gluon_frag->Draw();
    Quark_frag->Draw("SAME");
    legend_frag.Draw("SAME");
        
    c1->SaveAs(q_title_frag);
    c1->Clear();

    //Jet Major axis
    char q_title_major_axis_CL[80] ="";
    sprintf(q_title_major_axis_CL,"Jet_major_axis.pdf(");
    double scale_q_CL_major = 1./Quark_major_axis_CL->Integral();
    double scale_g_CL_major = 1./Gluon_major_axis_CL->Integral();
        
    Quark_major_axis_CL->Scale(scale_q_CL_major);
    Gluon_major_axis_CL->Scale(scale_g_CL_major);
    Gluon_major_axis_CL->SetMaximum(2.5*Quark_major_axis_CL->GetMaximum());
    Gluon_major_axis_CL->Draw();
    Quark_major_axis_CL->Draw("SAME");
    legend_major_axis_CL.Draw("SAME");
        
    c1->SaveAs(q_title_major_axis_CL);
    c1->Clear();
    
    char q_title_major_axis[80] ="";
    sprintf(q_title_major_axis,"Jet_major_axis.pdf)");
    double scale_q_part_major_axis = 1./Quark_major_axis->Integral();
    double scale_g_part_major_axis = 1./Gluon_major_axis->Integral();
        
    Quark_major_axis->Scale(scale_q_part_major_axis);
    Gluon_major_axis->Scale(scale_g_part_major_axis);
    Gluon_major_axis->SetMaximum(2.5*Quark_major_axis->GetMaximum());
    Gluon_major_axis->Draw();
    Quark_major_axis->Draw("SAME");
    legend_major_axis.Draw("SAME");
        
    c1->SaveAs(q_title_major_axis);
    c1->Clear();

    //Jet Minor axis
    char q_title_minor_axis_CL[80] ="";
    sprintf(q_title_minor_axis_CL,"Jet_minor_axis.pdf(");
    double scale_q_CL_minor = 1./Quark_minor_axis_CL->Integral();
    double scale_g_CL_minor = 1./Gluon_minor_axis_CL->Integral();
        
    Quark_minor_axis_CL->Scale(scale_q_CL_minor);
    Gluon_minor_axis_CL->Scale(scale_g_CL_minor);
    Gluon_minor_axis_CL->SetMaximum(2.5*Quark_minor_axis_CL->GetMaximum());
    Gluon_minor_axis_CL->Draw();
    Quark_minor_axis_CL->Draw("SAME");
    legend_minor_axis_CL.Draw("SAME");
        
    c1->SaveAs(q_title_minor_axis_CL);
    c1->Clear();
    
    char q_title_minor_axis[80] ="";
    sprintf(q_title_minor_axis,"Jet_minor_axis.pdf)");
    double scale_q_part_minor_axis = 1./Quark_minor_axis->Integral();
    double scale_g_part_minor_axis = 1./Gluon_minor_axis->Integral();
        
    Quark_minor_axis->Scale(scale_q_part_minor_axis);
    Gluon_minor_axis->Scale(scale_g_part_minor_axis);
    Gluon_minor_axis->SetMaximum(2.5*Quark_minor_axis->GetMaximum());
    Gluon_minor_axis->Draw();
    Quark_minor_axis->Draw("SAME");
    legend_minor_axis.Draw("SAME");
        
    c1->SaveAs(q_title_minor_axis);
    c1->Clear();
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
