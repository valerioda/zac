#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

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
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TApplication.h>
#include <TLegend.h>
#include "TCanvas.h"
#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "TStyle.h"
#include "TChain.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TMarker.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooUniform.h"
#include "RooPolynomial.h"
#include "RooPlot.h"

using namespace RooFit;

using namespace std;

int getOptions(int argc, char** argv);
void usage();
int readConfig(char*);
const std::vector<vector<double>> ReadTier2(char*, int);
double calibratePeak(double, double );
pair<double, double> approximatePeakPosition(double, double, int);

TH1D* h_ene;
TH1D* hTotalResiduals;
vector<double> data;
int gfitType = 0;
int fitType = 0;
double glob_interv = 0.;
// If you add more fit types, remember to update also the description in `usage'
const char* fitTypeNames[6] = { "Three Functions + Constant", "Three Functions + Line",
				"Two Functions + Constant", "Two Functions + Line",
				"Gaussian + Constant", "Gaussian + Line"};
const char* fitFCNNames[4] = { "Chi Square minimization", "Binned Maximum Likelihood", "Unbinned Maximum Likelihood", "NULL"};
const char* fitStatus[14] = {"CONVERGED", "ignored", "command unreadable", "unknown command", "NOT CONVERGED",
			     "command is a request to read PARAMETER definitions", "SET INPUT command", "SET TITLE command",
			     "SET COVAR command", "reserved", "END command", "EXIT or STOP command", "RETURN command", "STARTING STATE"};

double gCalib = 1;
const int nChannels = 65;
int actualChannels = 0;
int filterType = 1; // 0 = Gauss, 1 = ZAC, 2 = FADC
int TP = 0;
//double rebinFactor = 2.5;
double rebinFactor = 4.0;

int channelsToRun[nChannels] = {0};
double maxFADC[nChannels] = {0.0};
string channelsMapping[nChannels] = {"NULL"};
int detectorType[nChannels] = {0};

const int nMaxPeaks = 7;
int nPeaks = nMaxPeaks;
double peaksEnergy[nMaxPeaks] = {2614.511, 2103.512, 1620.738, 1592.513, 860.564, 727.33, 583.187};
double peaksEnergyErr[nMaxPeaks] = {0.010, 0.010, 0.010, 0.010, 0.03, 0.01, 0.01};
int peaksProcessed[nMaxPeaks] = {0,0,0,0,0,0,0};
int nChs_zac[nMaxPeaks]    = {150, 100, 50, 100, 80, 80, 80}; // zac
int nChs_gauss[nMaxPeaks]  = {100,  60, 50,  50, 40, 40, 40}; // gauss
int mapPeakstoFitType[nMaxPeaks] = {0, 3, 4, 2, 5, 5, 1};
double histo_norm = 1;
bool tier2_filename_as_is = false;

int fcnType = 2; // 1 = binned fit, 2 = unbinned fit
string tier2_dir = "/nfs/gerda5/gerda-data/blind/v06.02/gen/tier2/ged/cal/run0095/";
string tier2_name = "";
string configName = "";
char *resultDir;
int runNumber = -1;
vector<double> filterParameters; // = {Sigma, FlatTop, Tau, Length}
bool doAll = true;

int readConfig(char* configName) {//, vector<char*> channelsMapping, vector<double> maxFADC, vector<int> channelsToRun, vector<int> detectorType) {

  ifstream cfg;
  cfg.open(configName);
  string str1;
  int maxCh = 0;
  while(!cfg.eof()){ // To get you all the lines.
    getline(cfg,str1); // Saves the line in STRING.
    if(!cfg.good()) break;
    size_t stop1 = str1.find_first_of(" ");
    string str2 = str1.substr(stop1+1);
    int detNum = atoi(str1.substr(0,stop1).c_str());
    
    size_t stop2 = str2.find_first_of(" ");
    string str3 = str2.substr(stop2+1);
    string detName = str2.substr(0, stop2);
    
    size_t stop3 = str3.find_first_of(" ");
    string str4 = str3.substr(stop3+1);
    int detType = atoi(str3.substr(0, stop3).c_str());
    
    int chStat = atoi(str4.c_str());
    
    channelsMapping[detNum] = detName;
    channelsToRun[detNum] = chStat;
    detectorType[detNum] = detType;
    if(detNum > maxCh) maxCh = detNum;
  }
  cfg.close();
  
  return maxCh+1;
}

int main(int argc, char *argv[]){
  
  if(argc < 3) {
    usage();
    return 0;
  }
  
  if (getOptions(argc, argv)){
    cout << "ERROR: could not understand the arguments\n" << endl;
    usage();
    return 1;
  }
  
  //if (TP) rebinFactor = 20;
  const int maxCh = readConfig((char*)configName.c_str());
  cout << " Using tier2 keylist: " << tier2_name << endl;
  cout << " Save results in the directory: " << resultDir << endl;
  cout << endl << " Reading configuration from file " << configName << endl;
  cout << " Run n." << runNumber  << endl;
  cout << " Energy filter is: ";
  if (filterType==1) cout << "ZAC" << endl;
  else if (filterType==2) cout << "DPLMS" << endl;
  else cout << "GAUSS" << endl;
  for(unsigned int index=0; index<maxCh; ++index) {
    cout << "ch n." << index << " - " << channelsMapping[index];
    if(detectorType[index] == 0 ) cout << " (Enriched BEGe)";
    else if(detectorType[index] == 1 ) cout << " (Enriched Coaxial)";
    else if(detectorType[index] == 2 ) cout << " (Inverted Coaxial)";
    if(channelsToRun[index]) cout << " status: ON " << endl;
    else  cout << " status: OFF " << endl;
  }
  cout << endl;
  cout << " Total channels to run: " << maxCh << endl;
  
  cout << "Fit " << fitTypeNames[gfitType] << " using " << fitFCNNames[fcnType] << endl;
  
  int peaksToRun = nPeaks;
  if(!doAll && !TP) {
    peaksToRun = 1;
    cout << "Fit will be done only for the 2614 keV peak! " << endl;
  }
  else if (TP){
    peaksToRun = 1;
    cout << "Fit will be done only for the pulser! " << endl;
  }
  
  int index = 0;
  
  vector<double> residuals;
  
  ofstream json_lin;
  ofstream json_quad;

  if (doAll) {
    json_lin.open("linearCalibration.json", ofstream::out);
    json_quad.open("quadraticCalibration.json", ofstream::out);
    if( filterType > 0 ) {
      json_lin << "{\"ZAC\": {" << std::flush;
      json_quad << "{\"ZAC\": {" << std::flush;
    }
    //    else if(filterType==2) {
    // json_lin << "{\"DPLMS\": {" << std::flush;
    //  json_quad << "{\"DPLMS\": {" << std::flush;
    //}
    else {
      json_lin << "{\"Gauss\": {" << std::flush;
      json_quad << "{\"Gauss\": {" << std::flush;
    }
  }
  vector<vector<double>> energyVec = ReadTier2((char*)tier2_name.c_str(), runNumber);
  if (energyVec.size() == 0)
    return 1;

  for(int j=0; j<maxCh; j++) {

    if( filterType != 0 )
      maxFADC[j] = 50000;
    else
      maxFADC[j] = 20000;
    
    if (runNumber == 119) {
      maxFADC[j] = 10000;
      rebinFactor = 1.0;
    }
    if(channelsToRun[j] == 0) {
      cout << endl <<  "******************* channel n." << j << " is OFF *******************" << endl << endl;
      char thisLinCalib0[200];
      sprintf(thisLinCalib0, "\"%d\": { \"func\": \"[0]+[1]*x\", \"params\": {\"0\": %2.9e, \"1\": %2.9e}}", j, 0, 0);
      if(doAll) {
	json_lin << thisLinCalib0 << std::flush;
	if(j != maxCh-1) json_lin << ", " << std::flush;
      }
      continue;
    }
    
    int vSize = energyVec[j].size();
    cout << "events in Spectrum for channel n." << j << " : " << vSize << endl;
    
    if(vSize < 100) {
      cout << endl <<  "******************* Skipping channel n." << j << " *******************" << endl << endl;
      char thisLinCalib0[200];
      sprintf(thisLinCalib0, "\"%d\": { \"func\": \"[0]+[1]*x\", \"params\": {\"0\": %2.9e, \"1\": %2.9e}}", j, 0, 0);
      if(doAll) {
	json_lin << thisLinCalib0 << std::flush;
	if(j != maxCh-1) json_lin << ", " << std::flush;
      }
      continue;
    }

    TMarker* m[nMaxPeaks];
    cout << "******************* Analizing channel n." << j << " *******************" << endl;

    char sch[5];
    if(j<10) sprintf(sch,"ch");
    else sprintf(sch,"ch");
    
    char* name_stream;
    /*
    sprintf(name_stream, "./par_linear_%s%i.txt", sch, j);
    ofstream ofs_lin;
    ofs_lin.open(name_stream, ofstream::out);

    sprintf(name_stream, "./par_quadratic_%s%i.txt", sch, j);
    ofstream ofs_quad;
    ofs_quad.open(name_stream, ofstream::out);

    sprintf(name_stream, "./res_linear_%s%i.txt", sch, j);
    ofstream res_lin;
    res_lin.open(name_stream, ofstream::out);

    sprintf(name_stream, "./res_quadratic_%s%i.txt", sch, j);
    ofstream res_quad;
    res_quad.open(name_stream, ofstream::out);

    sprintf(name_stream, "./ene_linear_%s%i.txt", sch, j);
    ofstream ene_lin;
    ene_lin.open(name_stream, ofstream::out);

    sprintf(name_stream, "./ene_quadratic_%s%i.txt", sch, j);
    ofstream ene_quad;
    ene_quad.open(name_stream, ofstream::out);
    */
    if( filterType ==1)
      name_stream = Form("%s/ZAC-FWHM_chn%i.txt", resultDir, j);
    else 
      name_stream = Form("%s/FWHM_chn%i.txt", resultDir, j);
    ofstream fwhm_file;
    fwhm_file.open(name_stream, ofstream::out);
    
    char hStitle[20];
    sprintf(hStitle,"Spectrum Ch %i - %s", j, channelsMapping[j].c_str());
    char hSname[20];
    sprintf(hSname,"hSpectrum_ch%i", j);
    
    TH1D* hSpectrum;
    if( filterType != 0 )
      hSpectrum = new TH1D(hSname, hStitle, (maxFADC[j]-1)/rebinFactor, 1, maxFADC[j]);
    else hSpectrum = new TH1D(hSname, hStitle, maxFADC[j]-1, 1, maxFADC[j]);
    for (unsigned int i=0; i<energyVec[j].size(); i++) {
      //if ( j==9 ) cout << energyVec[j].at(i) << endl;
      double thisVal = energyVec[j].at(i);
      hSpectrum->Fill(thisVal);
    }
    //if ( j == 9 ){
    //cout << "BIN CONTENT CHN " << j << endl; 
    //for ( int i = 0; i < maxFADC[j]; i++ )
    //cout << hSpectrum->GetBinContent(i+1) << endl;
    //}

    double xEne[peaksToRun];
    double sqrtXEne[peaksToRun];
    double xEne_err[peaksToRun];

    double yFWHM[peaksToRun];
    double yFWHM_err[peaksToRun];

    double yPeak[peaksToRun];
    double yPeak_err[peaksToRun];

    Int_t  x, y;
    //UInt_t w=800;
    //UInt_t  height=600;
    //gVirtualX->GetWindowSize(gClient->GetRoot()->GetId(),x,y,w,height);
    TCanvas *c;
    char canvasName[20];
    sprintf(canvasName, "c_ch%i", j);
    
    if(peaksToRun>9) {
      c = new TCanvas(canvasName, canvasName, 800, 600);
      c->Divide(4,3);
    }
    else {
      c = new TCanvas(canvasName, canvasName, 600, 600);
      if ( doAll ) c->Divide(3,3);
    }
    for ( int np = 0; np < peaksToRun; np++ ) {
      TPad* thisPad; 
      c->cd(np+1);
      sprintf(canvasName, "c_ch%i_%i", j, np);
      char canvasTitle[30];
      sprintf(canvasTitle,"energy - ch.%i @ %2.2f",j, peaksEnergy[np]);
      if (doAll){
	thisPad = (TPad*)c->GetPad(np+1);
	thisPad->cd();
	thisPad->SetTitle(canvasTitle);
	thisPad->SetName(canvasName);
      }
      //else c->cd();
      if(np==0) {
	//if(detectorType[j] == 0) fitType = 0;
	//else fitType = 1;
	if (TP) fitType = 4;
	else fitType = gfitType;
      }
      else {
	fitType = mapPeakstoFitType[np];
      }
      cout << endl;
      cout << "#############################################" << endl;
      cout << "#                                           #" << endl;
      if(!TP) {
	cout << "#   Fitting peak at: " << peaksEnergy[np] << endl;
      }
      else {
	cout << "#   Fitting Pulser "<< endl;
      }
      cout << "#   with: " << fitTypeNames[fitType] << endl;
      cout << "#                                           #" << endl;
      cout << "#############################################" << endl;
      c->cd(np+1);
      gPad->SetGrid();
      pair<double,double> searchLimits;

      if(np==0) {
	sort(energyVec[j].begin(), energyVec[j].end());
	searchLimits.second = min(energyVec[j].at(vSize-1), maxFADC[j]); // Never go further than maxFADC[j]
	if (TP) searchLimits.first = energyVec[j].at(0);
	else searchLimits.first = searchLimits.second/3;
	cout << "first = " << searchLimits.first << endl;
	cout << "second = " << searchLimits.second << endl;
      }
      else searchLimits = approximatePeakPosition(peaksEnergy[np], gCalib, np);

      TH1D *h = new TH1D("h","", 300, searchLimits.first, searchLimits.second);
      h->Sumw2();
      for (unsigned int i=0; i<vSize; i++) {
	double thisEne = energyVec[j].at(i);
	if(thisEne <= searchLimits.second && thisEne >= searchLimits.first ) {
	  h->Fill(thisEne);
	  data.push_back(thisEne);
	}
      }
      double posPeak = searchLimits.second;
      //for (int i = 0; i < h->GetNbinsX(); i++) {
      //if( h->GetBinContent(h->GetNbinsX()-i) > 100 ) {
      //posPeak = 0.95*h->GetBinCenter(h->GetNbinsX()-i);
      //break;
      //}
      //}
      cout << "channel " << j << " max of the spectrum " << posPeak << endl;
      
      //--------- find parameter of the peak from spectrum -------
      h->GetXaxis()->SetRangeUser(searchLimits.first,posPeak);
      double average = h->GetBinCenter( h->GetMaximumBin());
      cout << "mean = " << average << endl;
      //double min = average - nChs[np];//* 0.99;
      //double max = average + nChs[np];//* 1.01;
      if (np==0) glob_interv = average*0.012;
      //else glob_interv *= 3.;
      double min = average - glob_interv;
      double max = average + glob_interv;
      if (TP){
	min = average - 10000;
	max = average + 10000;
      }
      h->Delete();

      int minIndex = -1;
      int maxIndex = 0;
      sort(data.begin(), data.end());
      for(int jj=0; jj<data.size(); jj++) {
	if(data[jj]>=min && minIndex<0) {
	  minIndex = jj;
	}
	if(data[jj]>max) {
	  maxIndex = data.size()-1-jj;
	  break;
	}
      }
      data.erase(data.begin(), data.begin() + minIndex);
      data.erase(data.end() - maxIndex, data.end());

      min = data.front();
      max = data.back();

      double maxH = max-min;
      int bin = (int)maxH;
      if( filterType != 0 ) {
	double temp = (double)bin/rebinFactor;//histo_norm;
	bin = (int)temp;
      }

      char henetitle[70];
      sprintf(henetitle,"Ch %i - %s - Peak @ %2.1f keV", j, channelsMapping[j].c_str(), peaksEnergy[np]);
      char henename[30];
      sprintf(henename,"h_ene_ch%i_E%2.0f",j,peaksEnergy[np]);

      h_ene = new TH1D(henename, henetitle, bin, 0, (int)maxH);
      cout << "nbin: " << bin << endl << "min = " << min << " max = " << max << " -> (0," << max-min << ")" << endl;
      if ( bin == 0 ) continue;
      peaksProcessed[np] = 1;
      h_ene->Sumw2();
      for (unsigned int i=0; i<energyVec[j].size(); i++) {
	double thisVal = energyVec[j].at(i) - min;
	if(thisVal >0 && thisVal<(maxH))
	  h_ene->Fill(thisVal);
      }
      for (unsigned int i=0; i<data.size(); i++)
	data[i] -= min;
      c->cd(np+1);
      gPad->SetLogy();
      c->Update();
      
      
      double A = h_ene->GetMaximum();

      // Initial value of the mean.  For the 2.6 MeV peak, look for the global
      // maximum of the full histogram.  In other cases look for the maximum in
      // a small region around peaksEnergy[np]/gCalib.  This avoids issues with
      // two close peaks in the same histogram.
      double thismean = 0.0;
      if (np==0)
	thismean = h_ene->GetBinCenter(h_ene->GetMaximumBin());
      else
	{
	  double approx_mean = peaksEnergy[np]/gCalib-min;
	  double range = ( filterType == 0 ) ? 10 : 20;
	  h_ene->GetXaxis()->SetRange(approx_mean-range, approx_mean+range);
	  thismean = h_ene->GetBinCenter(h_ene->GetMaximumBin());
	  h_ene->GetXaxis()->SetRange(0, max-min);
	}

      double sigma = 0.001 * (thismean + min);
      cout << "preliminary maximum = " << thismean << " A = " << A << endl;

      RooRealVar energy("x", "Energy", data.front(), data.back());
      RooRealVar mean("mean", "Mean", thismean, thismean - 2.0 *sigma, thismean + 2.0 * sigma);
      RooRealVar sgm("sgm", "Sigma", sigma, 0.3 * sigma, 10.0 * sigma);
      RooRealVar tailSigma("tailSigma", "Tail Sigma", sigma, 0.7 * sigma, 5.0 * sigma);
      // The peaks at 583, 727 and 860 keV are along a descent, so the slope should be negative.
      RooRealVar bkgSlope("bkgSlope", "bkgSlope", (np<3) ? 1e-3 : -1e-3, -0.01, 0.01);
      RooFormulaVar xx("xx", "x-mean", RooArgList(energy, mean));

      RooGaussian gauss("gauss","gaussian PDF", energy, mean, sgm);
      RooGenericPdf complErrorFun("complErrorFun", "complErrorFun",
				  "0.5 * TMath::Erfc((x - mean) / sgm / TMath::Sqrt(2.0))",
				  RooArgSet(energy, mean, sgm));
      RooGenericPdf errorExpFun("errorExpFun", "errorExpFun",
				"0.5 * TMath::Exp((x - mean) / tailSigma + sgm * sgm / (2 * tailSigma * tailSigma)) * TMath::Erfc((x - mean) / sgm / TMath::Sqrt(2.0) + sgm / tailSigma / TMath::Sqrt(2.0))",
				RooArgSet(energy, mean, sgm, tailSigma));
      RooUniform constantBkg("constbkg", "constbkg", RooArgSet());
      RooPolynomial linearBkg("linearbkg", "linearbkg", xx, RooArgList(bkgSlope));

      if (TP)
	{
	  sgm = 10.0 * sigma;
	  sgm.setRange(0.0, 10000.*sigma);
	}

      double flatbkg = 1;
      int nn = 0;
      for (int kk=3*bin/4; kk<=bin; kk++) {
	flatbkg += h_ene->GetBinContent(kk);
	nn++;
      }
      flatbkg *= (max-min)/nn/4.0; // The average is computed in the last quarter

      RooRealVar n1("n1", "Component #1 fraction", data.size() / 2, data.size() / 10, 2 * data.size());
      RooRealVar n2("n2", "Component #2 fraction", 1, 0, 2 * data.size());
      RooRealVar n3("n3", "Component #3 fraction", 1, 0, data.size() / 2);
      RooRealVar n4("n4", "Component #4 fraction", flatbkg, 0, flatbkg * 6);
      RooRealVar n5("n5", "Component #5 fraction", flatbkg, 0, flatbkg * 6);
      RooExtendPdf e1("e1", "extended pdf #1", gauss, n1);
      RooExtendPdf e2("e2", "extended pdf #2", complErrorFun, n2);
      RooExtendPdf e3("e3", "extended pdf #3", errorExpFun, n3);
      RooExtendPdf e4("e4", "extended pdf #4", constantBkg, n4);
      RooExtendPdf e5("e5", "extended pdf #5", linearBkg, n5);

      vector<RooRealVar> parameters;

      RooPlot* xframe = energy.frame();
      xframe->SetTitle(henetitle);
      //sgm.getVal()
      //int rangeu = (int)mean.getVal()-min;
      //xframe->GetXaxis()->SetRangeUser(rangeu-10,rangeu+10);
      xframe->GetXaxis()->SetTitle(Form("FADC-%1.2f",min));
      xframe->GetYaxis()->SetTitle("counts");
      xframe->GetXaxis()->CenterTitle();
      xframe->GetYaxis()->CenterTitle();

      RooAddPdf *model;

      if (fitType == 0)
	model = new RooAddPdf("model", fitTypeNames[fitType], RooArgList(e1, e2, e3, e4));
      else if(fitType == 1)
	model = new RooAddPdf("model", fitTypeNames[fitType], RooArgList(e1, e2, e3, e5));
      else if(fitType == 2)
	model = new RooAddPdf("model", fitTypeNames[fitType], RooArgList(e1, e2, e4));
      else if(fitType == 3)
	model = new RooAddPdf("model", fitTypeNames[fitType], RooArgList(e1, e2, e5));
      else if(fitType == 4)
	model = new RooAddPdf("model", fitTypeNames[fitType], RooArgList(e1, e4));
      else if(fitType == 5)
	model = new RooAddPdf("model", fitTypeNames[fitType], RooArgList(e1, e5));

      cout << "INITIAL VALUES: Mean = " << mean.getVal() << " Sigma = " << sgm.getVal() <<  " Tail Sigma = " << tailSigma.getVal() << " Bkg = ";
      if ((fitType%2) == 0)
	cout << n4.getVal() << endl;
      else
	cout << n5.getVal() << endl;

      if(fcnType==1)
	{
	  RooDataHist channel_energy("channel_energy", "channel_energy", energy, Import(*h_ene)) ;
	  cout << " ******************************* USING RooFit FOR THE BINNED FIT" << endl;
	  model->fitTo(channel_energy);
	  //channel_energy.plotOn(xframe, MarkerColor(kRed+1), LineColor(kAzure+2), MarkerStyle(8), LineWidth(1));
	  channel_energy.plotOn(xframe, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(7), LineWidth(2));
	}
      else if(fcnType==2)
	{
	  RooDataSet channel_energy("channel_energy", "channel_energy", RooArgSet(energy));
	  for (unsigned int i=0; i<data.size(); i++)
	    {
	      energy = data[i];
	      channel_energy.add(RooArgSet(energy));
	    }
	  cout << " ******************************* USING RooFit FOR THE UNBINNED FIT" << endl;
	  model->fitTo(channel_energy);
	  //channel_energy.plotOn(xframe, MarkerColor(kRed+1), LineColor(kAzure+2), MarkerStyle(8));
	  channel_energy.plotOn(xframe, MarkerColor(kBlue), LineColor(kBlue), MarkerStyle(7), LineWidth(2));
	}
      else
	{
	  cout << " ******************************* FITTING METHOD NOT SUPPORTED, ABORTING" << endl;
	  return 1;
	}
      //model->plotOn(xframe, LineColor(kSpring+3), LineStyle(1), LineWidth(2));
      model->plotOn(xframe, LineColor(kRed), LineStyle(1), LineWidth(3));

      parameters.push_back(mean);
      parameters.push_back(sgm);
      if (fitType == 0) parameters.push_back(tailSigma);
      if ((fitType%2) == 1) parameters.push_back(bkgSlope);
      cout << "FINAL VALUES: Mean = " << mean.getVal() << " Sigma = " << sgm.getVal() <<  " Tail Sigma = " << tailSigma.getVal() << " Bkg = ";
      if ((fitType%2) == 0)
	cout << n4.getVal() << endl;
      else 
	cout << n5.getVal() << endl;
      c->cd(np+1);
      xframe->Draw();
      m[np] = new TMarker(mean.getVal() + min, A, 4);
      c->Update();
      gStyle->SetOptStat(0);
      sqrtXEne[np] = sqrt(peaksEnergy[np]);
      xEne[np] = peaksEnergy[np];
      xEne_err[np] = peaksEnergyErr[np];
      double muCh = mean.getVal() + min;
      double muCh_err = mean.getError();
      yPeak[np] = muCh;
      yPeak_err[np] = muCh_err;
      double calib = calibratePeak(muCh,peaksEnergy[np]);
      if(np==0) gCalib = calib;
      if (TP) calib = 1.;
      double mu = muCh*calib;
      double mu_err = muCh_err*calib;
      double sigma_2 = sgm.getVal() * calib;
      double sigma_err = sgm.getError() * calib;
      double FWHM = 2.35*sigma_2;
      if (TP) FWHM = 2.35*sigma_2/mu;
      double FWHM_err = 2.35*sigma_err;
      if (TP) FWHM_err = 2.35*(sigma_err/sigma_2+mu_err/mu)*FWHM;
      yFWHM[np] = FWHM;
      yFWHM_err[np] = FWHM_err;
      double chi2Fit = xframe->chiSquare(parameters.size());
      //if(!doAll) {
      //char out_file[100];
      //sprintf(out_file,"./FWHM_chn%i.txt",j);
      //std::ofstream myFile;
      //myFile.open (out_file, std::ios_base::app|std::ios_base::out);
      //cout << parameters.size() << endl;
      //cout << filterParameters.size() << endl;
      if (isnormal(chi2Fit) && FWHM_err < 1){
	if (filterParameters.size()>0){
	  if (filterType==2)
	    fwhm_file << peaksEnergy[np] << " " << filterParameters.at(0) << " " << filterParameters.at(1) << " " << filterParameters.at(2) <<  " " << filterParameters.at(3) << " " << filterParameters.at(4) << " " << filterParameters.at(5) << " " << filterParameters.at(6) << " " << FWHM << " " << FWHM_err << " " << endl;//chi2Fit << " " << parameters.size() << endl;
	  else
	    fwhm_file << filterParameters.at(0) << " " << filterParameters.at(1) << " " << filterParameters.at(2) <<  " " << filterParameters.at(3) << " " << FWHM << " " << FWHM_err << " " << chi2Fit << " " << parameters.size() << endl;
	  //myFile.close();
	}
	else fwhm_file << j << " " << peaksEnergy[np] <<  " " << FWHM << " " << FWHM_err << " " << chi2Fit << " " << parameters.size() << endl;
      }
      else {
	if (filterParameters.size()>0){
	  if (filterType==2)
	    fwhm_file << peaksEnergy[np] << " 0 0 0 0 0 0 0 0 0" << endl;
	  else
	    fwhm_file << "0 0 0 0 0 0 0 0" << endl;
	}
	else fwhm_file << j << " " << peaksEnergy[np] <<  " 0 0 0 0 " << endl;
      }
      
      //-------------- print the result --------------------
      TPaveText *label = new TPaveText(0.6, 0.25, 0.9, 0.45,"NDC");

      char result[100];
      sprintf(result,"FWHM = %2.3f #pm %2.3f keV",FWHM,FWHM_err);
      label->AddText(result);
      for (int npar = 0; npar < parameters.size(); npar++) {
	if(npar==0)
	  sprintf(result,"%s = %2.3f #pm %2.3f", parameters[npar].GetName(), parameters[npar].getVal() + min, parameters[npar].getError());
	else
	  sprintf(result,"%s = %2.3f #pm %2.3f", parameters[npar].GetName(), parameters[npar].getVal(), parameters[npar].getError());
	label->AddText(result);
      }
      sprintf(result,"Chi2/ndof = %2.3f", chi2Fit);
      label->AddText(result);
      label->SetTextAlign(12);
      label->SetTextSize(0.025);
      label->SetBorderSize(1);
      label->SetFillColor(0);
      c->cd(np+1);
      label->Draw();
      //h_ene->GetXaxis()->SetRangeUser(mean.getVal()-3*sgm.getVal(), mean.getVal()+3*sgm.getVal());
      c->Update();

      data.clear();
      h_ene->Reset();
      delete model;
    } //loop on peaks

    c->cd();
    char *fileoutname;
    if(filterType == 1)
      fileoutname = Form("%s/ZAC-Fit_ch%i.pdf", resultDir, j);
    else
      fileoutname = Form("%s/Fit_ch%i.pdf", resultDir, j);
    c->Print(fileoutname);


    if(doAll) {
      TCanvas* cv = new TCanvas(1);
      cv->cd();
      TGraphErrors* gr_FWHMvsE = new TGraphErrors(peaksToRun, xEne, yFWHM, xEne_err, yFWHM_err);
      char grTitle[30];
      sprintf(grTitle,"Ch %i - %s", j, channelsMapping[j].c_str());
      gr_FWHMvsE->SetTitle(grTitle);
      gr_FWHMvsE->GetXaxis()->SetTitle("energy [keV]");
      gr_FWHMvsE->GetYaxis()->SetTitle("FWHM[keV]");
      gr_FWHMvsE->SetMarkerColor(kRed+1);
      gr_FWHMvsE->SetMarkerStyle(8);
      gr_FWHMvsE->GetXaxis()->SetRangeUser(-10,3000);
      gr_FWHMvsE->Draw("AP");
      cv->Update();

      int actualP = 0;
      double xxFWHM[nPeaks], yyFWHM[nPeaks], erryyFWHM[nPeaks];

      for (int nPoints = 0; nPoints<nPeaks; nPoints++ ) {
	if(nPoints == 1)
	  {
	    cout <<  "... skipping peak @ " << xEne[nPoints] << " from FWHM fit "<< endl;
	  }
	else
	  {
	    xxFWHM[actualP] = xEne[nPoints];
	    yyFWHM[actualP] = yFWHM[nPoints];
	    erryyFWHM[actualP] = yFWHM_err[nPoints];
	    actualP++;
	  }
      }

      TGraphErrors FWHMCalib = TGraphErrors(actualP, xxFWHM, yyFWHM, NULL, erryyFWHM);
      TF1 fwhmFitFunc = TF1("fwhmFitFunc","sqrt([0] + [1]*x)", -10, 5000);
      fwhmFitFunc.SetParNames("A", "B"); // The names for the parameters are needed when summing several TF1s with parameters
      fwhmFitFunc.SetParameters(6, 2e-3); // Initialize parameters to some reasonable values, needed for the fit
      FWHMCalib.Fit(&fwhmFitFunc);
      fwhmFitFunc.SetLineColor(kSpring+3);
      fwhmFitFunc.Draw("same");
      TF1 fwhmFitFuncError = TF1("fwhmFitFuncError","sqrt(([2] * 0.5 / sqrt([0] + [1]*x)) ** 2.0 + ([3] * 0.5 * x / sqrt([0] + [1]*x)) ** 2.0)", -10, 5000);
      fwhmFitFuncError.SetParNames("A", "B", "A_err", "B_err");
      fwhmFitFuncError.SetParameters(fwhmFitFunc.GetParameter(0), fwhmFitFunc.GetParameter(1),
				     fwhmFitFunc.GetParError(0), fwhmFitFunc.GetParError(1));
      TF1 fwhmFitFuncPlus = TF1("fwhmFitFuncPlus","fwhmFitFunc + fwhmFitFuncError", -10, 5000);
      fwhmFitFuncPlus.SetLineColor(kSpring+3);
      fwhmFitFuncPlus.SetLineStyle(2);
      fwhmFitFuncPlus.Draw("same");
      TF1 fwhmFitFuncMinus = TF1("fwhmFitFuncMinus","fwhmFitFunc - fwhmFitFuncError", -10, 5000);
      fwhmFitFuncMinus.SetLineColor(kSpring+3);
      fwhmFitFuncMinus.SetLineStyle(2);
      fwhmFitFuncMinus.Draw("same");
      cv->Update();

      //-------------- print the result --------------------
      TPaveText *fwhmlabel = new TPaveText(0.2, 0.5, 0.4, 0.8,"NDC");

      char fwhmresult[100];
      fwhmlabel->AddText("FWHM = #sqrt{A + B #upoint E[keV]}");
      sprintf(fwhmresult,"A = %2.2f #pm %2.2f", fwhmFitFunc.GetParameter(0), fwhmFitFunc.GetParError(0));
      fwhmlabel->AddText(fwhmresult);
      sprintf(fwhmresult,"B = %2.2f #pm %2.2f", fwhmFitFunc.GetParameter(1), fwhmFitFunc.GetParError(1));
      fwhmlabel->AddText(fwhmresult);
      sprintf(fwhmresult,"FWHM(0) = %2.2f #pm %2.2f", fwhmFitFunc.Eval(0), fwhmFitFuncError.Eval(0));
      fwhmlabel->AddText(fwhmresult);
      sprintf(fwhmresult,"FWHM(Q_{#beta#beta}) = %2.2f #pm %2.2f", fwhmFitFunc.Eval(2039), fwhmFitFuncError.Eval(2039));
      fwhmlabel->AddText(fwhmresult);

      sprintf(fwhmresult,"Chi2/ndof = %2.2f", FWHMCalib.Chisquare(&fwhmFitFunc) / (actualP - 2));
      fwhmlabel->AddText(fwhmresult);
      fwhmlabel->SetTextAlign(12);
      fwhmlabel->SetTextSize(0.025);
      fwhmlabel->SetBorderSize(0);
      fwhmlabel->SetFillColor(0);
      fwhmlabel->Draw();
      cv->Update();
      cv->Print(Form("%s/Fit_FWHM_chn%i.pdf", resultDir, j));

      cout << "Fit FWHM = sqrt((" << fwhmFitFunc.GetParameter(0)  << " +/- " << fwhmFitFunc.GetParError(0) <<
	") + x * (" << fwhmFitFunc.GetParameter(1) << " +/- " << fwhmFitFunc.GetParError(1) << "))" << endl;

      TCanvas* cv2 = new TCanvas(1);
      cv2->cd();
      TGraphErrors* gr_Calib = new TGraphErrors(peaksToRun, yPeak, xEne, yPeak_err, xEne_err);
      gr_Calib->SetTitle(grTitle);
      gr_Calib->GetYaxis()->SetTitle("energy [keV]");
      gr_Calib->GetXaxis()->SetTitle("#ch[a.u]");

      gr_Calib->SetMarkerColor(kSpring+3);
      gr_Calib->SetMarkerStyle(8);
      gr_Calib->Draw("AP");
      TF1* pol1 = new TF1("pol1","pol1",0, 3000);
      gr_Calib->Fit(pol1,"F");

      vector<double> xx, yy, erryy;
      for (int nPoints = 0; nPoints<nPeaks; nPoints++ ) {
	xx.push_back(xEne[nPoints]);
	yy.push_back(yPeak[nPoints]);
	erryy.push_back(yPeak_err[nPoints]);
      }

      TGraphErrors* gr_CalibTEMP = new TGraphErrors(peaksToRun, xEne, yPeak, xEne_err, yPeak_err);
      TF1* pol1_1;
      //if( zac == 1 ) pol1_1 = new TF1("pol1_1","[0]*x",0, 3000);
      //else
      pol1_1 = new TF1("pol1_1","pol1",0, 3000);
      gr_Calib->Fit(pol1_1,"");

      //***** uncomment for NON analytic fit *****
      //double Afit = (zac == 0) ? pol1_1->GetParameter(0) : 0;
      //double Bfit = (zac == 1) ? pol1_1->GetParameter(0) : pol1_1->GetParameter(1);
      //double errAfit = (zac == 0) ? pol1_1->GetParError(0) : 0;
      //double errBfit = (zac == 1) ? pol1_1->GetParError(0) : pol1_1->GetParError(1);
      //***************************************

      double Afit = pol1_1->GetParameter(0);
      double Bfit = pol1_1->GetParameter(1);
      double errAfit = pol1_1->GetParError(0);
      double errBfit = pol1_1->GetParError(1);

      cout << "Energy linear Fit -> A = " << Afit << " +/- " << errAfit << " B = " << Bfit << " +/- " << errBfit << endl;

      double eneChi2 = 0;
      //const int SS = yy.size();
      const int SS = peaksToRun;
      double thisRes[SS];
      double thisRes_err[SS];
      for (int ii=0; ii<SS; ii++) {
	double valF = Afit + Bfit*yPeak[ii];
	double valF_err = sqrt(errAfit*errAfit + yPeak[ii]*yPeak[ii]*errBfit*errBfit + Bfit*Bfit*yPeak_err[ii]*yPeak_err[ii]);
	//ene_lin << xEne[ii] << " " << valF << " " << valF_err << " " << yPeak[ii] << " " << yPeak_err[ii] << endl;
	//double val = yy[ii] - valF;
	double val = xEne[ii] - valF;
       	residuals.push_back(val);
	thisRes[ii] = val;
	thisRes_err[ii] =  sqrt(xEne_err[ii]*xEne_err[ii] + errAfit*errAfit + yPeak[ii]*yPeak[ii]*errBfit*errBfit + Bfit*Bfit*yPeak_err[ii]*yPeak_err[ii]);//erryy[ii];
	double thisChi2 = val*val;
	eneChi2 += thisChi2;
	//res_lin << xEne[ii] << " " << thisRes[ii] << " " << thisRes_err[ii] << endl;
	// //cout << " ----> " << ii << " " << xx[ii] << " "<< yy[ii] << " " << valF << " " << val << " " << eneChi2 << endl;
      }
      eneChi2 /= ((double)SS-2.);
      //eneChi2 =  pol1->GetChisquare()/pol1->GetNDF();

      TF1* eneFitFunc = new TF1("eneFitFunc","[0] + [1]*x",0,15000);
      eneFitFunc->SetParameter(0, Afit);
      eneFitFunc->SetParameter(1, Bfit);
      eneFitFunc->SetLineColor(kRed+1);
      eneFitFunc->Draw("same");
      TF1* eneFitFuncPlus = new TF1("eneFitFuncPlus","[0] + [1]*x + sqrt([2]*[2] + x*x*[3]*[3])",0,15000);
      eneFitFuncPlus->SetParameter(0, Afit);
      eneFitFuncPlus->SetParameter(1, Bfit);
      eneFitFuncPlus->SetParameter(2, errAfit);
      eneFitFuncPlus->SetParameter(3, errBfit);
      eneFitFuncPlus->SetLineColor(kRed+1);
      eneFitFuncPlus->SetLineStyle(2);
      eneFitFuncPlus->Draw("same");
      TF1* eneFitFuncMinus = new TF1("eneFitFuncMinus","[0] + [1]*x - sqrt([2]*[2] + x*x*[3]*[3])",0,15000);
      eneFitFuncMinus->SetParameter(0, Afit);
      eneFitFuncMinus->SetParameter(1, Bfit);
      eneFitFuncMinus->SetParameter(2, errAfit);
      eneFitFuncMinus->SetParameter(3, errBfit);
      eneFitFuncMinus->SetLineColor(kRed+1);
      eneFitFuncMinus->SetLineStyle(2);
      eneFitFuncMinus->Draw("same");
      cv2->Update();

      //-------------- print the result --------------------
      TPaveText *enelabel = new TPaveText(0.2, 0.5, 0.5, 0.8,"NDC");

      char eneresult[100];
      enelabel->AddText(" y = A + B #upoint x");
      sprintf(eneresult,"A = %2.2e #pm %2.2e",Afit, errAfit);
      enelabel->AddText(eneresult);
      sprintf(eneresult,"B = %2.2e #pm %2.2e",Bfit, errBfit);
      enelabel->AddText(eneresult);

      char thisLinCalib[200];
      sprintf(thisLinCalib, "\"%d\": { \"func\": \"[0]+[1]*x\", \"params\": {\"0\": %2.9e, \"1\": %2.9e}}", j, Afit, Bfit);
      if(doAll) {
	json_lin << thisLinCalib << std::flush;
	if(j != maxCh-1) json_lin << ", " << std::flush;
      }

      //ofs_lin << Afit << " " << errAfit << " " << Bfit << " " << errBfit  << endl;
      sprintf(eneresult,"Chi2/ndof = %2.3f", eneChi2);
      enelabel->AddText(eneresult);
      enelabel->SetTextAlign(12);
      enelabel->SetTextSize(0.025);
      enelabel->SetBorderSize(0);
      enelabel->SetFillColor(0);
      enelabel->Draw();
      cv2->Update();
      cv2->Print(Form("%s/Fit_Energy_linear_ch%i.pdf", resultDir, j));

      TCanvas* cv3 = new TCanvas(1);
      cv3->cd();
      TGraphErrors* gr_Calib2 = new TGraphErrors(peaksToRun, yPeak, xEne, yPeak_err, xEne_err);
      gr_Calib2->SetTitle(grTitle);
      gr_Calib2->GetYaxis()->SetTitle("energy [keV]");
      gr_Calib2->GetXaxis()->SetTitle("#ch[a.u]");
      gr_Calib2->SetMarkerColor(kAzure-2);
      gr_Calib2->SetMarkerStyle(8);
      gr_Calib2->Draw("AP");
      TF1* pol2;
      pol2 = new TF1("pol2","pol2",0, 3000);
      gr_Calib2->Fit(pol2,"");


      double p0 = 0;
      double p1 = 0;
      double p2 = 0;
      double p0err = 0;
      double p1err = 0;
      double p2err = 0;

      p0 = pol2->GetParameter(0);
      p1 = pol2->GetParameter(1);
      p2 = pol2->GetParameter(2);
      p0err = pol2->GetParError(0);
      p1err = pol2->GetParError(1);
      p2err = pol2->GetParError(2);
      TPaveText *enelabel2 = new TPaveText(0.2, 0.5, 0.5, 0.8,"NDC");

      char eneresult2[100];
      enelabel2->AddText(" y = A + B #upoint x + C #upoint x^{2}");
      sprintf(eneresult2,"A = %2.2e #pm %2.2e",p0, p0err);
      enelabel2->AddText(eneresult2);
      sprintf(eneresult2,"B = %2.2e #pm %2.2e",p1, p1err);
      enelabel2->AddText(eneresult2);
      sprintf(eneresult2,"C = %2.2e #pm %2.2e",p2, p2err);
      enelabel2->AddText(eneresult2);

      sprintf(eneresult2,"Chi2/ndof = %2.3f/%d = %2.3f", pol2->GetChisquare(), pol2->GetNDF(), pol2->GetChisquare()/pol2->GetNDF() );
      enelabel2->AddText(eneresult2);
      enelabel2->SetTextAlign(12);
      enelabel2->SetTextSize(0.025);
      enelabel2->SetBorderSize(0);
      enelabel2->SetFillColor(0);
      enelabel2->Draw();
      cv3->Update();

      char thisQuadCalib[200];
      sprintf(thisQuadCalib, "\"%d\": { \"func\": \"[0]+[1]*x+[2]*x*x\", \"params\": {\"0\": %2.9e, \"1\": %2.9e, \"2\": %2.9e}}", j, p0, p1, p2);
      if(doAll) {
	json_quad << thisQuadCalib << std::flush;
	if(j != maxCh-1) json_quad << ", " << std::flush;
      }

      //ofs_quad << p0 << " " << p0err << " " << p1 << " " << p1err << " " << p2 << " " << p2err << endl;

      cv3->Print(Form("%s/Fit_Energy_quadratic_ch%i.pdf", resultDir, j));
      double thisRes2[SS];
      double thisRes2_err[SS];
      for (int ii=0; ii<yy.size(); ii++) {
	double valF = p0 + p1*yPeak[ii] + p2*yPeak[ii]*yPeak[ii];
	double valF_err = p0err*p0err + p1err*p1err*yPeak[ii]*yPeak[ii] + yPeak[ii]*yPeak[ii]*yPeak[ii]*yPeak[ii]*p2err*p2err + (p1 + 2.0*yPeak[ii]*p2)*(p1 + 2.0*yPeak[ii]*p2)*yPeak_err[ii]*yPeak_err[ii];
	double val = xEne[ii] - valF;
	thisRes2[ii] = val;

	thisRes2_err[ii] = sqrt(xEne_err[ii]*xEne_err[ii]
				+ (p1 + 2.0*yPeak[ii]*p2)*(p1 + 2.0*yPeak[ii]*p2)*yPeak_err[ii]*yPeak_err[ii]);
	//ene_quad << xEne[ii] << " "<< valF << " " << valF_err << endl;
      	//res_quad << xEne[ii] << " "<< thisRes2[ii] << " " << thisRes2_err[ii] << endl;
      }
      TCanvas* cSpectrum = new TCanvas(2);
      cSpectrum->cd();
      cSpectrum->SetLogy();
      hSpectrum->GetXaxis()->SetTitle("#FADC");
      hSpectrum->GetYaxis()->SetTitle("counts");
      hSpectrum->Draw();
      for(int mm=0; mm<nPeaks; mm++) {
	if (!peaksProcessed[mm]) continue;
	m[mm]->SetMarkerColor(kRed+1);
	m[mm]->Draw();
	char strLab[20];
	sprintf(strLab,"%2.1f keV",peaksEnergy[mm]);
	TLatex *energyLabel = new TLatex(0.9*m[mm]->GetX(), 1.2*m[mm]->GetY(),strLab);
	energyLabel->SetTextFont(42);
	energyLabel->SetTextColor(kRed+1);
	energyLabel->SetTextSize(0.03);
	energyLabel->Draw("same");
      }
      cSpectrum->Update();
      cSpectrum->Print(Form("%s/Spectrum_ch%i.pdf", resultDir,j));
      TCanvas* cResiduals = new TCanvas(1);
      cResiduals->cd();

      TGraphErrors* gr_Residuals = new TGraphErrors(SS, xEne, thisRes, 0, thisRes_err);
      gr_Residuals->SetTitle(grTitle);
      gr_Residuals->GetYaxis()->SetTitle("residuals [keV]");
      gr_Residuals->GetYaxis()->SetRangeUser(-1, 1);
      gr_Residuals->GetXaxis()->SetTitle("energy [keV]");
      gr_Residuals->SetMarkerColor(kRed-4);
      gr_Residuals->SetMarkerStyle(8);
      gr_Residuals->Draw("AP");

      TGraphErrors* gr_Residuals2 = new TGraphErrors(SS, xEne, thisRes2, 0, thisRes2_err);
      gr_Residuals2->SetTitle(grTitle);
      gr_Residuals2->SetMarkerColor(kAzure);
      gr_Residuals2->SetMarkerStyle(4);
      gr_Residuals2->Draw("P");

      TLegend* leg = new TLegend(0.25,0.75,0.5,0.85);
      leg->SetLineColor(0);
      leg->AddEntry(gr_Residuals,"Linear Fit","pl");
      leg->AddEntry(gr_Residuals2,"Quadratic Fit","pl");
      leg->Draw();

      cResiduals->Update();
      cResiduals->Print(Form("%s/Residuals_ch%i.pdf", resultDir, j));
      
      index++;
      //delete hSpectrum;
    } //if(doAll)

    //ofs_lin.close();
    //ofs_quad.close();
    //res_lin.close();
    //res_quad.close();
    //ene_lin.close();
    //ene_quad.close();
    fwhm_file.close();

    cout << "*****************************************************************************" << endl;
  }

  if(doAll) {
    json_lin << "}}" << std::flush;
    json_lin.close();
    json_quad << "}}" << std::flush;
    json_quad.close();
  }
  return 0;
}

//______________________________________________________________________________

int getOptions(int argc, char **argv){
  int c;
  while ((c = getopt(argc, argv, "c:D:d:f:F:m:o:r:t:T:z:")) != -1)
    switch (c) {
    case 'c':
      configName = string(optarg);
      break;
    case 'D':
      resultDir = optarg;
      break;
    case 'd':
      tier2_dir = string(optarg);
      break;
    case 'f':
      gfitType = atoi(optarg);
      break;
    case 'F':
      filterType = atoi(optarg);
      break;
    case 'm':
      fcnType = atoi(optarg);
      break;
    case 'o':
      tier2_filename_as_is = true;
      break;
    case 'r':
      runNumber = atoi(optarg);
      break;
    case 't':
      tier2_name = string(optarg);
      break;
    case 'T':
      TP = atoi(optarg);
      break;
    case 'z':
      string delimiter = ",";
      string s = string(optarg) + delimiter;
      size_t pos = 0;
      for( int i = 0; i < 4; i++){
	pos = s.find(delimiter);
	if (pos != string::npos){
	  filterParameters.push_back(atof(s.substr(0, pos).c_str()));
	  s.erase(0, pos + delimiter.length());
	}
	else
	  return 1;
      }
      doAll = false;
      break;
    }
  if (TP)
    doAll = false;
  return 0;
}

//______________________________________________________________________________

void usage() {
  cout << "USAGE: ./energyCalibration [options]" << endl;
  cout << "List of options:" << endl;
  cout << " -c: configuration file" << endl;
  cout << " -D: result directory" << endl;
  cout << " -d: directory of the tier2 file" << endl;
  cout << " -t: tier2 file name" << endl;
  cout << " -f: fit type (";
  for ( int n = 0; n < 3; n++)
    {
      if (n > 0) cout << "               ";
      cout << 2*n   << " = " << fitTypeNames[2*n] << "; ";
      cout << 2*n+1 << " = " << fitTypeNames[2*n+1];
      if (n==2) cout << ")";
      cout << endl;
    }
  cout << " -F: filter type (0 = Gauss, 1 = ZAC (default), 2 = DPLMS)" << endl;
  cout << " -m: fit mode (1 = binned (default), 2 = unbinned)" << endl;
  cout << " -o: if this option is used, the files paths in the file list will be read as they are" << endl;
  cout << " -r: run number" << endl;
  cout << " -T: whether to fit pulser only (0 = do not fit pulser (default), 1 = fit pulser only)" << endl;
  cout << " -z: comma-separated list of filter parameters Sigma,FlatTop,Tau,Length" << endl;
}

//______________________________________________________________________________

const std::vector<vector<double>> ReadTier2(char* tier2_name, int run) {
  
  std::vector<vector<double>> energyV;
  ifstream file;
  file.open(tier2_name);
  
  int nFiles = 0;
  while (1){
    char file_name[200];
    file.getline(file_name, 200);
    if (!file.good()) break;
    TFile* tierfile = NULL;
    if(tier2_filename_as_is)
      tierfile = TFile::Open(Form("%s", file_name));
    else
      tierfile = TFile::Open(Form("%s-ged-tier2.root", tier2_dir.c_str(), file_name));
    TTree* tier2 = (TTree*)tierfile->Get("tier2");
    std::vector<double>* gaussEnergy = NULL;
    if(tier2->GetBranch("GEMDEnergyGauss_energy"))
      tier2->SetBranchAddress("GEMDEnergyGauss_energy", &gaussEnergy);
    std::vector<double>* energy = NULL;
    if(filterType == 1)
      //tier2->SetBranchAddress("GEMDGenericShaping_energy", &energy);
      tier2->SetBranchAddress("GEMDZACShaping_energy", &energy);
    else if(filterType == 2)
      //tier2->SetBranchAddress("GEMDFADC_energy", &energy);
      tier2->SetBranchAddress("GEMDDPLMSShaping_energy", &energy);
    else
      energy = gaussEnergy;
    std::vector<int>* waveformTag = NULL;
    bool waveformBranch = false;
    if(tier2->GetBranch("GEMDFADC_waveformTag")) {
      tier2->SetBranchAddress("GEMDFADC_waveformTag", &waveformTag);
      waveformBranch = true;
    }
    std::vector<int>* quality = NULL;
    bool qualityBranch = false;
    if(tier2->GetBranch("GEMDQuality_isProcessed")) {
      tier2->SetBranchAddress("GEMDQuality_isProcessed", &quality);
      qualityBranch = true;
    }
    std::vector<double>* fitExpCoefficient = NULL;
    bool fitExpCoefficientBranch = false;
    if(tier2->GetBranch("GEMDBaseline_fitExpCoefficient")) {
      tier2->SetBranchAddress("GEMDBaseline_fitExpCoefficient", &fitExpCoefficient);
      fitExpCoefficientBranch = true;
    }
    std::vector<double>* baselineSigma = NULL;
    bool baselineSigmaBranch = false;
    if(tier2->GetBranch("GEMDBaseline_baselineSigma")) {
      tier2->SetBranchAddress("GEMDBaseline_baselineSigma", &baselineSigma);
      baselineSigmaBranch = true;
    }

    std::vector<double>* auxBaseline1Sigma = NULL;
    bool auxBaseline1SigmaBranch = false;
    if(tier2->GetBranch("GEMDBaseline_auxBaseline1Sigma")) {
      tier2->SetBranchAddress("GEMDBaseline_auxBaseline1Sigma", &auxBaseline1Sigma);
      auxBaseline1SigmaBranch = true;
    }

    std::vector<double>* triggerTime = NULL;
    bool triggerBranch = false;
    if(tier2->GetBranch("GEMDTrigger_trigger")) {
      tier2->SetBranchAddress("GEMDTrigger_trigger", &triggerTime);
      triggerBranch = true;
    }
    std::vector<int>* triggerNumber = NULL;
    bool triggerNumberBranch = false;
    if(tier2->GetBranch("GEMDFTTrigger_triggerNumber")) {
      tier2->SetBranchAddress("GEMDFTTrigger_triggerNumber", &triggerNumber);
      triggerNumberBranch = true;
    }

    std::vector<double>*  risetime = NULL;
    bool risetimeBranch = false;
    if(tier2->GetBranch("GEMDRiseTimeHF_risetime")) {
      tier2->SetBranchAddress("GEMDRiseTimeHF_risetime", &risetime);
      risetimeBranch = true;
    }

    std::vector<double>*  maxAmpTime = NULL;
    bool maxAmpTimeBranch = false;
    if(tier2->GetBranch("GEMDEnergyGauss_maxAmpTime")) {
      tier2->SetBranchAddress("GEMDEnergyGauss_maxAmpTime", &maxAmpTime);
      maxAmpTimeBranch = true;
    }

    std::vector<double>*  energyRevPol = NULL;
    bool energyRevPolBranch = false;
    if(tier2->GetBranch("GEMDEnergyGauss_energyRevPol")) {
      tier2->SetBranchAddress("GEMDEnergyGauss_energyRevPol", &energyRevPol);
      energyRevPolBranch = true;
    }
    std::vector<double>* baseline = NULL;
    bool baselineBranch = false;
    if(tier2->GetBranch("GEMDBaseline_baseline")) {
      tier2->SetBranchAddress("GEMDBaseline_baseline", &baseline);
      baselineBranch = true;
    }
    std::vector<int>*  eventType = NULL;
    tier2->SetBranchAddress("GEMDFADC_eventType", &eventType);
    tier2->GetEntry(0);
    int channels = energy->size();
    UInt_t n_events = tier2->GetEntries();
    cout << "n. events = " << n_events << endl;
    if(nFiles==0) energyV.resize(channels);
    for(UInt_t i = 0; i < n_events; i++) {
      tier2->GetEntry(i);
      for(int channel = 0; channel < channels; channel++){
	// Use the `is0vbbFromCal' cuts in
	// https://github.com/mppmu/gerda-metadata/blob/master/config/_aux/cuts/default-cuts.json
	
	bool isWaveformGood = false;
	if(waveformBranch) {
	  if(waveformTag->at(channel) == 0)
	    isWaveformGood = true;
	}
	else
	  isWaveformGood = true;
	bool isQualityGood = false;
	if(qualityBranch) {
	  if(quality->at(channel) == 1)
	    isQualityGood = true;
	}
	else
	  isQualityGood = true;

	bool isBaselineSlopeGood = false;
	if(fitExpCoefficientBranch) {
	  if(fitExpCoefficient->size() == 0 ||
	     TMath::Abs(fitExpCoefficient->at(channel)) <= 200)
	    isBaselineSlopeGood = true;
	}
	else
	  isBaselineSlopeGood = true;
	
	bool isPrePileUpGood = false;
	if(baselineSigmaBranch && auxBaseline1SigmaBranch) {
	  if(baselineSigma->at(channel)/auxBaseline1Sigma->at(channel) < 1.3)
	    isPrePileUpGood = true;
	}
	else
	  isPrePileUpGood = true;
	bool isTriggerGood = false;
	if(triggerBranch) {
	  double thisTriggerTime = triggerTime->at(channel);
	  if(thisTriggerTime > 77e3 && thisTriggerTime < 84e3)
	    isTriggerGood = true;
	}
	else
	  isTriggerGood = true;
	bool isTriggerNumberGood = false;
	if(triggerNumberBranch) {
	  if(triggerNumber->at(channel) == 1)
	    isTriggerNumberGood = true;
	}
	else
	  isTriggerNumberGood = true;
	bool isRisetimeGood = false;
	if(baselineSigmaBranch && risetimeBranch) {
	  double thisRisetime = risetime->at(channel);
	  if((gaussEnergy->at(channel) < 3*baselineSigma->at(channel) ||
	      thisRisetime > 70.0) && thisRisetime < 5000.0)
	    isRisetimeGood = true;
	}
	else
	  isRisetimeGood = true;
	bool isMaxAmpTimeGood = false;
	if(baselineSigmaBranch && risetimeBranch) {
	  double thismaxAmpTime = maxAmpTime->at(channel);
	  if(thismaxAmpTime > 142e3 && thismaxAmpTime < 148e3)
	    isMaxAmpTimeGood = true;
	}
	else
	  isMaxAmpTimeGood = true;
	bool isNegPolarityGood = false;
	if(energyRevPolBranch && baselineSigmaBranch && baselineBranch) {
	  if(TMath::Abs(energyRevPol->at(channel)) <=
	     6*baselineSigma->at(channel) ||
	     baseline->at(channel) < -65535+100)
	    isNegPolarityGood = true;
	}
	else
	  isNegPolarityGood = true;
	bool complexCut = false;
	int thisEventType = eventType->at(channel);
	bool isGaussEnergyGood = false;
	if(tier2->GetBranch("GEMDEnergyGauss_energy"))
	  isGaussEnergyGood = gaussEnergy->at(channel) > 1;
	if(TP)
	  {
	    if(thisEventType)
	      complexCut = true;
	  }
	else if(!thisEventType && isWaveformGood && isQualityGood )//&&
	  //isBaselineSlopeGood && isPrePileUpGood && isTriggerGood &&
	  //isTriggerNumberGood )//&& isRisetimeGood && isMaxAmpTimeGood &&
	  //isGaussEnergyGood && isNegPolarityGood)
	  complexCut = true;
	if(complexCut) {
	  if(i==0)
	    actualChannels = channels;
	  energyV[channel].push_back(energy->at(channel));
	}
      }
    }
    tierfile->Close();
    delete tierfile, tier2, gaussEnergy, energy, waveformTag, quality,
      fitExpCoefficient, baselineSigma, auxBaseline1Sigma, triggerTime, risetime,
      maxAmpTime, energyRevPol, baseline, eventType;
    nFiles++;
  }
  file.close();
  if (nFiles == 0)
    {
      cout << "ERROR: Could not fine any key in " << tier2_name << endl;
      return energyV;
    }

  return energyV;
}

double calibratePeak( double fadc, double energyPeak)
{
  return energyPeak/fadc;
}

//given a 2614 calibration returns the boudary search for others peaks
pair<double,double> approximatePeakPosition(double peakenergy, double calibration, int np)
{
  pair<double,double> searchLimits;
  double fadc = peakenergy/calibration;
  double range;
  if (filterType==0)
    range = nChs_gauss[np];
  else
    range = nChs_zac[np];
  // We need to enlarge the range on the left to make it symmetric.
  // TODO: discover why this is necessary.
  searchLimits.first = fadc - range * 1.5;
  searchLimits.second = fadc + range;
  cout << peakenergy << "/" << calibration << " -> "<< searchLimits.first << " - " << searchLimits.second << endl;
  return searchLimits;
}
