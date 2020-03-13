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

	// Lorenzian Peak function
Double_t dtermFunc(Double_t *x, Double_t *par) {
double M_p=0.938;
    return -(2.768)/( (1+(x[0]/(1.204*1.204)) )*(1+(x[0]/(1.204*1.204))) );
	
}

/*Double_t dtermFunc1(Double_t *x, Double_t *par) {
double M_p=0.938;
    return (-2)/( (1+(x[0]/(0.84*0.84)) )*(1+(x[0]/(0.84*0.84))) );
	
}
*/

Double_t dtermFunc1(Double_t *x, Double_t *par) {
double M_p=0.938;
    return  -2*(1./((1 + x[0]/0.71)*(1 + x[0]/0.71)))*(1/(1 + x[0]/(4.*M_p*M_p) ))*( 1 + 2.79*x[0]/( 4*M_p*M_p ) );;
	
}


int DtermPlot() {
gROOT->SetBatch(kTRUE);
	

 auto *g = new TGraph("gpd_dterm.dat");
  auto *g1 = new TGraph("gpd_dterm1.dat");

 double x_old[17]={0.0000,0.050,0.1000,0.1500,0.200,0.2500,0.3000,0.3500,0.4000,0.4500,0.5000,0.5500,0.6000,0.6500,0.7000,0.7500,0.8000};
    double y_old[17]={-3.28889,-2.86802, -2.52306, -2.23679,-1.99661,-1.79314,-1.61926,-1.46949,-1.33958 ,-1.22617,-1.12657,-1.03864,-0.96061,-0.89105,-0.82879,-0.77283,-0.72235};

int n = g->GetN(); double *x = g->GetX(); double *y = g->GetY(); 
for (int i=0;i<n;i++) {  y[i] = y[i]*(1-x[i]); }


  auto *dtermOld = new TGraph(17,x_old,y_old);

TCanvas *canvas  = new TCanvas("canvas");
//gPad->SetLogy();



g->Draw("AL");
//g1->Draw("L same");
cout<<(g->Integral())<<endl;

canvas->SaveAs("comparesimudata.pdf"); 

double integ=(g->Integral());   

TF1 *Dterm = new TF1("Dterm",dtermFunc,0.0,0.8,0);
TF1 *Dterm1 = new TF1("Dterm",dtermFunc1,0.0,0.8,0);

TCanvas *canvas1  = new TCanvas("canvas1");
//gPad->SetLogy();




dtermOld->Draw();
Dterm->Draw("same");
Dterm1->SetLineColor(kBlue);
Dterm1->Draw("same");
//g1->Draw("L same");
// cout<<(g1->Integral())<<endl;

canvas1->SaveAs("dterm.pdf"); 


	
return 123456789;
}
