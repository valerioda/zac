#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include <TApplication.h>
#include <TLegend.h>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TPaveLabel.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TMarker.h"

using namespace std;

int main(int argc, char *argv[]){
  //gROOT->Reset();
  
  //  TApplication *myapp = new TApplication("myapp",0,0);
  
  char* scanDir = argv[1];//scan path
  int run = atoi(argv[2]);//number of run
  const int nChn = atoi(argv[3]);//number of total channels
  char* date = argv[4];//calibration date
  char* time = argv[5];
  
  int firstChn = 24;
  
  int ThreeDPlot = 0;
  if (argc > 6)
    ThreeDPlot = atoi(argv[6]);
  int save_pdf = 1;
  if (argc > 7)
    save_pdf = atoi(argv[7]);
  double plot = 0.;
  if (argc > 8)
    plot = atof(argv[8]);
  double chi_lim = 100.;
  if (argc > 9)
    chi_lim = atof(argv[9]);
  
  //-------------------------------------------------------//
  
  char JSON_file[200];
  if ( run < 10) sprintf(JSON_file,"%s/gerda-run000%d-%sT%sZ-cal-ged-tier2-calib.json",scanDir,run,date,time);
  else if ( run < 100) sprintf(JSON_file,"%s/gerda-run00%d-%sT%sZ-cal-ged-tier2-calib.json",scanDir,run,date,time);
  else if ( run < 1000) sprintf(JSON_file,"%s/gerda-run0%d-%sT%sZ-cal-ged-tier2-calib.json",scanDir,run,date,time);
  else sprintf(JSON_file,"%s/gerda-run%d-%sT%sZ-cal-ged-tier2-calib.json",scanDir,run,date,time);
  std::ofstream myJSON;
  myJSON.open(JSON_file);
  myJSON << "{" << endl;
  myJSON << "  \"gelatio\": { " << endl;
  myJSON << "    \"env\": {" << endl;
      
  char out_file[200];
  sprintf(out_file,"%s/ZACanalysis.txt",scanDir);
  std::ofstream myFile;
  myFile.open(out_file);
  
  char ini_file[200];
  sprintf(ini_file,"%s/optZAC_%s.ini",scanDir,date);
  std::ofstream myFile2;
  myFile2.open(ini_file);
  
  
  //--------------------------------------------------//
  TFile *grafico = new TFile(Form("%s/grafico_%s.root",scanDir,date),"recreate") ;

  double minFWHM[nChn];
  double minFWHM_err[nChn];
  double minSigma[nChn];
  double minTau[nChn];
  double minFT[nChn];

  for (int chn = firstChn; chn < nChn; chn++){
    char *name = Form("%s/ZAC-FWHM_channel%d.txt",scanDir,chn);
    //double ft_val = 1.;
    //double tau_val = 160.;
    ifstream file;
    file.open(name);
    if(file) {
      double x1,x2,x3,x4;
      double x5,x6,x7,x8;
      double sigma[500];
      double tau[500];
      double FT[500];
      double FT_1[500];
      double FWHM[500];
      double FWHM_err[500];
      double chi2[500];
      double sigma_1[500];
      double FWHM_1[500];
      double FWHM_1_err[500];
      double tau_1[500];
      int n_lines = 0, n_lines_1 = 0, j = 0;
      while (1){
	file >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8;
	if (!file.good()) break;
	sigma[n_lines] = x1;
	FT[n_lines] = x2;
	tau[n_lines] = x3;
	FWHM[n_lines] = x5;
	FWHM_err[n_lines] = x6;
	chi2[n_lines] = fabs(x7/x8);
	n_lines++;
      }
      file.close();
      
      minFWHM[chn] = FWHM[0];
      minFWHM_err[chn] = FWHM_err[0];
      minSigma[chn] = 10;
      minTau[chn] = 160;
      minFT[chn] = 1;
      TH2D *sigmaft = new TH2D("sigmaft","", 6, 0, 30, 5, 0.5, 3);
      for (int i = 0; i < n_lines; i++){
	if ((chi2[i] < chi_lim)){// && (tau[i]==tau_val) && (FWHM_err[i]<0.06)){
	  FT_1[j] = FT[i];
	  sigma_1[j] = sigma[i];
	  FWHM_1[j] = FWHM[i];
	  FWHM_1_err[j] = FWHM_err[i];
	  tau_1[j] = tau[i];
	  int ind1 = sigma_1[j]/5;
	  int ind2 = FT_1[j]*2-1;
	  sigmaft->SetBinContent( ind1, ind2, FWHM_1[j]);
	  n_lines_1++;
	  if (FWHM_1[j] < minFWHM[chn]){
	    minFWHM[chn] = FWHM_1[j];
	    minFWHM_err[chn] = FWHM_1_err[j];
	    minSigma[chn] = sigma_1[j];
	    minTau[chn] = tau_1[j];
	    minFT[chn] = FT_1[j];
	  }
	  j++;
	}
      }
      
      
      /////////////////////////////////////////////////
      TCanvas *c1 = new TCanvas(Form("c%d",chn),"c1");
      
      if ( n_lines_1 > 0 ) {
	//c1->SetGrid();
	//TGraph2D *g_1_3d;
	TGraph *g_1_2d;
	if ( ThreeDPlot ) {
	  //g_1_3d = new TGraph2D( j, FT_1, sigma_1, FWHM_1 );
	  //g_1_3d->Draw("COLZ2");
	  sigmaft->SetTitle("");
	  sigmaft->GetXaxis()->SetTitle("Sigma [#mus]");
	  sigmaft->GetYaxis()->SetTitle("Flat Top [#mus]");
	  //sigmaft->GetZaxis()->SetTitle("FWHM [keV]");
	  //sigmaft->GetXaxis()->SetTitleSize(0.045);
	  //sigmaft->GetYaxis()->SetTitleSize(0.045);
	  //sigmaft->GetXaxis()->SetNdivisions(3);
	  //sigmaft->GetYaxis()->SetNdivisions(4);
	  sigmaft->Draw("COLZ2");
	  sigmaft->SetStats(0);
	  
	  //g_1_3d->SetTitle(Form(" %s - ZAC optimization - chn %d", date, chn));
	  //g_1_3d->GetXaxis()->SetTitle("Sigma [#mus]");
	  //g_1_3d->GetYaxis()->SetTitle("FT [#mus]");
	  //g_1_3d->GetZaxis()->SetTitle("FWHM [keV]");
	  //g_1_3d->GetXaxis()->CenterTitle();
	  //g_1_3d->GetYaxis()->CenterTitle();
	  //g_1_3d->GetZaxis()->CenterTitle();
	  //g_1_3d->GetZaxis()->SetRange(2, 4);
	  c1->Update();
	  if (save_pdf) c1->Print(Form("%s/ch%d_FWHMvsSigma-FT.pdf",scanDir,chn-firstChn));
	}
	else if( plot == 1 ) {
	  g_1_2d = new TGraphErrors (j,FT_1,FWHM_1,0,FWHM_1_err);
	  g_1_2d->SetMarkerColor(4);
	  g_1_2d->SetLineColor(4);
	  g_1_2d->SetTitle(Form(" %s - ZAC optimization - chn %d", date, chn-firstChn));
	  g_1_2d->GetXaxis()->SetTitle("FT [#mus]");
	  g_1_2d->GetYaxis()->SetTitle("FWHM [keV]");
	  if (chn < 10) g_1_2d->SetName(Form("g_2d_chn0%d", chn-firstChn));
	  else g_1_2d->SetName(Form("g_2d_chn%d", chn-firstChn));
	  //g_1_2d->GetXaxis()->CenterTitle();
	  //g_1_2d->GetYaxis()->CenterTitle();
	  //	g_1_2d->GetYaxis()->SetRangeUser(2.4, 6);
	  g_1_2d->Draw("APE");
	  //TFile *grafico = TFile::Open(Form("grafico_%s.root",date));
	  grafico->cd();
	  g_1_2d->Write();
	  //grafico->Close();
	  c1->Update();
	  if (save_pdf) c1->Print(Form("%s/ch%d_FWHMvsFT.pdf",scanDir,chn-firstChn));
	}
	else if( plot == 2 ) {
	  g_1_2d = new TGraphErrors (nChn, tau, FWHM, 0,FWHM_err);
	  g_1_2d->SetMarkerColor(4);
	  g_1_2d->SetLineColor(4);
	  g_1_2d->SetTitle(Form(" %s - ZAC optimization - chn %d", date, chn-firstChn));
	  g_1_2d->GetXaxis()->SetTitle("tau [#mus]");
	  g_1_2d->GetYaxis()->SetTitle("FWHM [keV]");
	  if (chn < 10) g_1_2d->SetName(Form("g_2d_chn0%d", chn-firstChn));
	  else g_1_2d->SetName(Form("g_2d_chn%d", chn-firstChn));
	  g_1_2d->GetXaxis()->CenterTitle();
	  g_1_2d->GetYaxis()->CenterTitle();
	  //	g_1_2d->GetYaxis()->SetRangeUser(2.4, 6);
	  g_1_2d->Draw("APE");
	  //TFile *grafico = TFile::Open(Form("grafico_%s.root",date));
	  grafico->cd();	
	  g_1_2d->Write();
	  //grafico->Close();
	  c1->Update();
	  if (save_pdf) c1->Print(Form("%s/ch%d_FWHMvsTau.pdf",scanDir,chn-firstChn));
	}
	else {
	  g_1_2d = new TGraphErrors (j,sigma_1,FWHM_1,0,FWHM_1_err);
	  g_1_2d->SetMarkerColor(4);
	  g_1_2d->SetLineColor(4);
	  g_1_2d->SetTitle(Form(" %s - ZAC optimization - chn %d", date, chn-firstChn));
	  g_1_2d->GetXaxis()->SetTitle("Sigma [#mus]");
	  g_1_2d->GetYaxis()->SetTitle("FWHM [keV]");
	  if (chn < 10) g_1_2d->SetName(Form("g_2d_chn0%d", chn-firstChn));
	  else g_1_2d->SetName(Form("g_2d_chn%d", chn-firstChn));
	  //g_1_2d->GetXaxis()->CenterTitle();
	  //g_1_2d->GetYaxis()->CenterTitle();
	  //	g_1_2d->GetYaxis()->SetRangeUser(2.4, 6);
	  g_1_2d->Draw("APE");
	  //grafico = TFile::Open(Form("grafico_%s.root",date));
	  grafico->cd();	
	  g_1_2d->Write();
	  c1->Update();
	  if (save_pdf) c1->Print(Form("%s/ch%d_FWHMvsSigma.pdf",scanDir,chn-firstChn));  
	}
      }
    }
    else {
      minFWHM[chn] = -1;
      minFWHM_err[chn] = -1;
      minSigma[chn] = 10;
      minTau[chn] = 160;
      minFT[chn] = 1;
    }
    
  }
  grafico->Close();
  
  for (int chn = firstChn; chn < nChn; chn++){
    char out[200];
    char ini[200];
    char ini2[200];
    char JSON_out[200];
    char JSON_out2[200];
    
    
    if (chn < 10)
      sprintf(out,"%d  %4.2f %4.2f %g %3.0f %.0f \n",chn-firstChn, minFWHM[chn], minFWHM_err[chn], minSigma[chn], minFT[chn], minTau[chn]);
    else
      sprintf(out,"%d %4.2f %4.2f %g %3.0f %.0f \n",chn-firstChn, minFWHM[chn], minFWHM_err[chn], minSigma[chn], minFT[chn], minTau[chn]);
    sprintf(ini,"FilterName[%d]=/nfs/gerda2/LNGSMiB/ZAC_filters/ZACfilter_L155_sigma%.0f_FT%g_tau%.0f.txt",chn-firstChn, minSigma[chn], minFT[chn], minTau[chn]);
    sprintf(ini2,"PulserFilterName[%d]=/nfs/gerda2/LNGSMiB/ZAC_filters/ZACfilter_L155_sigma%.0f_FT%g_tau%.0f.txt",chn-firstChn, minSigma[chn], minFT[chn], minTau[chn]);
    //if ( chn >= 36 && run > 93) {
    if ( chn >= 60 && run > 93) {
      if (chn == (nChn-1))
	sprintf(JSON_out,"      \"zac_filter_physig_%d\": \"L155_sigma%.0f_FT1p5_tau%.0f\"",chn-firstChn, minSigma[chn], minTau[chn]);
      else
	sprintf(JSON_out,"      \"zac_filter_physig_%d\": \"L155_sigma%.0f_FT1p5_tau%.0f\",",chn-firstChn, minSigma[chn], minTau[chn]);
      
      if (chn == (nChn-1))
	sprintf(JSON_out2,"      \"zac_filter_pulser_%d\": \"L155_sigma%.0f_FT1p5_tau%.0f\"",chn-firstChn, minSigma[chn], minTau[chn]);
      else
	sprintf(JSON_out2,"      \"zac_filter_pulser_%d\": \"L155_sigma%.0f_FT1p5_tau%.0f\",",chn-firstChn, minSigma[chn], minTau[chn]);
    } 
    else {
      if (chn == (nChn-1))
	sprintf(JSON_out,"      \"zac_filter_physig_%d\": \"L155_sigma%.0f_FT1_tau%.0f\"",chn-firstChn, minSigma[chn], minTau[chn]);
      else
	sprintf(JSON_out,"      \"zac_filter_physig_%d\": \"L155_sigma%.0f_FT1_tau%.0f\",",chn-firstChn, minSigma[chn], minTau[chn]);
      
      if (chn == (nChn-1))
	sprintf(JSON_out2,"      \"zac_filter_pulser_%d\": \"L155_sigma%.0f_FT1_tau%.0f\"",chn-firstChn, minSigma[chn], minTau[chn]);
      else
	sprintf(JSON_out2,"      \"zac_filter_pulser_%d\": \"L155_sigma%.0f_FT1_tau%.0f\",",chn-firstChn, minSigma[chn], minTau[chn]);
    }
    
    myFile << out;
    
    myFile2 << ini << endl;
    //myFile2 << ini2 << endl;
    
    myJSON << JSON_out << endl;
    //myJSON << JSON_out2 << endl;
  }
  
  //myJSON.open(JSON_file,std::ios_base::app|std::ios_base::out);
  myJSON << "    }" << endl;
  myJSON << "  } " << endl;
  myJSON << "}" << endl;
  
  myJSON.close();
  myFile.close();
  myFile2.close();
  //myapp->Run();
  return 0;
}
  

