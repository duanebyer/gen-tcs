#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraphPainter.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TF1.h"
#include "TTreeReader.h"
#include "TEfficiency.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH3F.h"


using namespace std;


int CheckCFFPlots() {
gROOT->SetBatch(kTRUE);
	

//gStyle->SetOptStat(111);
	gStyle->SetPalette(55);	
	//gStyle->SetLabelSize(.04, "xyz");
	gStyle->SetTitleSize(.05, "xyz");
	gStyle->SetTitleSize(.07,"t");
	gStyle->SetLineWidth(1);
	gStyle->SetFrameLineWidth(1); 
	gStyle->SetGridStyle(3);
		gStyle->SetGridColor(kGray);
gStyle->SetGridWidth(1);
	/*
	gStyle->SetHistLineWidth(2);*/
	//gStyle->SetMarkerStyle(13);

TCanvas *canvas1 = new TCanvas("c11","multipads", 2000, 1400);
TCanvas *canvas = new TCanvas("c1","multipads", 2000, 1400);
  //gStyle->SetOptTitle(0);
    //gStyle->SetLabelOffset(0.2);
  canvas->Divide(6,6);


ifstream infile("CFF_Check.dat", ios::binary);
	


    

    Int_t  nv=0;
    Double_t e_vc=0.;
    Double_t Q2=0., t=0., xi=0., ImH=0.0 , ReH=0.0 , ImE=0.0, ReE=0.0, ImHt=0.0, ReHt=0.0, D=0.0;

    for(int i=0;i<6;i++) { //q2
        for(int j=0;j<6;j++) { // xi
        
       TPad *p1 = new TPad(Form("ampp%d",(i+100*j)), "p1", 0.05+0.15*i, 0.05+0.15*j, 0.2+0.15*i, 0.2+0.15*j, 0, 0, 0);
 p1->SetTopMargin(0);
    p1->SetRightMargin(0);
    if(j>0)p1->SetBottomMargin(0);
        if(i>0)p1->SetLeftMargin(0);
  p1->Draw();
     p1->SetGridy();
          p1->SetGridx();
  
        TGraphErrors *g1 = new TGraphErrors();
    g1->SetName(Form("amp%d",(i+100*j)));
    
    g1->SetLineColor(kRed);
        
        	for(int k=0;k<11;k++) { //t
        	
            if(!infile.good()) break;

            infile >> Q2 >> xi >> t >> ImH >> ReH >> ImE>> ReE>> ImHt >>ReHt>> D;

            g1->SetPoint(k, -t, ImH);
            
            cout<<t<<" "<<ImH<<" "<<endl;
            g1->SetTitle(Form("xi %f Q2 %f  ;-t ;ImH",xi,Q2));
           
            }
           // g1->GetXaxis()->SetRangeUser(0.,0.6);
           // g1->GetHistogram()->GetXaxis()->Set(3,0.0,0.6);
            /*canvas->cd(1+i+6*j);*/p1->cd(); g1 -> Draw("AL");
            cout<<(1+i+6*j)<<endl;
             canvas->cd();
             //canvas1->cd();g1 -> Draw("AL");
           //g1->Delete();
        }
    }
    
    
    infile.close();
    
/* 
 
  
  p1->cd(); gr -> Draw("AP");
  p2->cd(); gr -> Draw("AP");*/



canvas->SaveAs("CheckCFF.pdf"); 
canvas1->SaveAs("CheckCFF1.pdf"); 

	
return 123456789;
}
