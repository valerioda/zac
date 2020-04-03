#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

void usage(){

  cout << "usage: ZACfilter <binning> <length> <cusp par> <flat top> <tau>" << endl;
}

int main( int argc, const char* argv[])
{
  
  //cout << argc << endl;

  if(argc<7) {
    usage();
    return 1;
  }
  const char *resDir = argv[1];
  const double bin = atoi(argv[2]);
  
  //parameters value in mu
  const double  N_t = atof(argv[3]); //filter lenght
  const double sigma_t = atof(argv[4]); //cusp par
  const double FT_t = atof(argv[5]); //flat top
  const double tau_t = atof(argv[6]); //preamplifier tau
  int daq = 0;
  if (argc>7) daq = atoi(argv[7]); //preamplifier tau
  
  if( daq == 1 ) cout << "FlashCam Data" << endl;
  else cout << "Struck Data" << endl;
  
  cout << "Binning of DAQ: " << bin << " ns" << endl;
  cout << "ZAC parameters in mus(N,sigma,FT,tau): " 
       << N_t << ", " << sigma_t << ", " << FT_t 
       << ", " << tau_t << endl;
  
  
  char *filename = Form("%s/ZACfilter_L%g_sigma%g_FT%g_tau%g.txt",resDir,N_t,sigma_t,FT_t,tau_t);
  if (argc>8) filename = Form("%s",argv[8]);
 
  ofstream file(filename); 
  ofstream file2(Form("%s/ZAC.txt",resDir)); 
  double FT_ns = FT_t*1000.;
  
  int N = N_t*1000./(double)bin;
  int sigma = sigma_t*1000./(double)bin;
  int FT = FT_ns/(double)bin;
  int tau = tau_t*1000./(double)bin;
  int L = 0;
  if ( (N-FT) % 2 == 0 )
    L = (N-FT)/2.;
  else
    {
      L = (N-FT-1)/2.;
      FT += 1;
    }

  ///////////////////////////////////////
  //creation of the cusp-like filter
  double filtro_cusp[N];
  double t[N];
  for (int i = 0; i < L; i++)
    {
      t[i] = i/1e3*bin;
      filtro_cusp[i] = sinh( ((double)i+1.)/(double)sigma) / sinh((double)L/(double)sigma );
      //filtro_cusp[N-1-i] = filtro_cusp[i];
    }
  for (int i = L; i < L+FT; i++)
    {
      t[i] = i/1e3*bin;
      filtro_cusp[i] = 1.;
    }
    for (int i = L+FT; i < N; i++)
      {
	t[i] = i/1e3*bin;
	int j = N - i - 1;
	filtro_cusp[i] = filtro_cusp[j];//sinh( ((double)j+1.)/(double)sigma) / sinh((double)L/(double)sigma );
      }
    
  //////////////////////////////////////////////
  //creation of the negative paraboles
  double P[N] ;
  for (int i = 0; i < L; i++)
    {
      P[i] = ((double)i + 1. -(double)L*0.5)*((double)i + 1. - (double)L*0.5) -((double)L*0.5)*((double)L*0.5);
      P[N-i] = P[i];
    }
  for (int i = L+FT; i < N; i++)
    {
      int j = N-i-1;
      P[i] = P[j];//((double)j + 1. -(double)L*0.5)*((double)j + 1. - (double)L*0.5) -((double)L*0.5)*((double)L*0.5);
    }
  
  double sum_P = 0;
  double sum_cusp = 0;
  for  (int i = 0; i < N; i++){
    sum_P = sum_P + P[i];
    sum_cusp = sum_cusp + filtro_cusp[i];}

  cout << "Area Cuspide = " << sum_cusp << endl;
  
  for  (int i = 0; i < N; i++)
    P[i] = -P[i]/sum_P*sum_cusp;//this is to normalize the area
  
  /////////////////////////////////////////////
  //creation of ZAC filter
  double zac[N];
  for  (int i = 0; i < N; i++)
    {
      zac[i]=filtro_cusp[i]+P[i];
    }

  //////////////////////////////////////
  //creation of deconvolution filter
  double filtro_dec[2];
  filtro_dec[0] = 1;
  filtro_dec[1] = -exp(-1./tau);
  //convolution process
  double tmp;
  double filtro[N];
  
  for (int i = 0; i < N; i++){
    int i1 = i;
    tmp = 0.;
    for (int j = 0; j < 2; j++){
      if(i1 >= 0 && i1 < N)
	tmp = tmp + (zac[i1]*filtro_dec[j]);
      i1 = i1-1;
    }
    //if ( daq == 1 ) filtro[i] = tmp;
    //else filtro[i] = -tmp;
    filtro[i] = tmp;
  }
  for (int i = 1; i < N; i++){
    file << filtro[i] << endl;
    file2 << zac[i] << endl;
  }
  file.close();
  file2.close();
  
  ///////////////////////////////////////////////
  //graphs of filter
  TGraph *g_cusp = new TGraph (N,t,filtro_cusp);
  TGraph *g_par = new TGraph (N,t,P);
  TGraph *g_zac = new TGraph (N,t,zac);
  g_zac->SetTitle("");
  g_zac->GetYaxis()->SetTitle("amplitude");
  g_zac->GetXaxis()->SetTitle("time [#mus]");
  //g_zac->GetXaxis()->CenterTitle();
  //g_zac->GetYaxis()->CenterTitle();
  g_zac->GetXaxis()->SetRangeUser(0,160);
  g_zac->SetLineWidth(2);
  g_zac->SetLineColor(2);
  g_cusp->SetLineWidth(2);
  g_cusp->SetLineStyle(2);
  g_cusp->SetLineColor(3);
  g_par->SetLineWidth(2);
  g_par->SetLineColor(4);
  g_par->SetLineStyle(2);
    
  TCanvas *c1 = new TCanvas("c1","ZAC filter");
  //c1->SetGrid();
  g_zac->Draw("AL");
  g_par->Draw("L");
  g_cusp->Draw("L");
  g_zac->Draw("L");
  TLegend *legend = new TLegend(.75,.80,.95,.95);
  legend->AddEntry(g_zac,"ZAC filter","l");
  legend->AddEntry(g_cusp,"cusp","l");
  legend->AddEntry(g_par,"parabolas","l");
  legend->Draw();
  c1->Print(Form("%s/ZAC.pdf",resDir));  
  
  //graph of deconvolved filter
  TGraph *g_dec = new TGraph (N,t,filtro);
  g_dec->SetTitle("");
  g_dec->GetYaxis()->SetTitle("amplitude");
  g_dec->GetXaxis()->SetTitle("time [#mus]");
  g_dec->GetXaxis()->SetRangeUser(0,160);
  g_dec->GetYaxis()->SetTitleOffset(1.5);
  g_dec->SetLineWidth(2);
  g_dec->SetLineColor(kBlue);
  //g_dec->SetMarkerSize(4);
  g_dec->SetMarkerColor(4);
  
  TCanvas *c2 = new TCanvas("c2","ZAC filter");
  //c2->SetGrid();
  c2->cd();
  //g_zac->Draw("AL");
  g_dec->Draw("AL");
  TLegend *legend1 = new TLegend(.75,.80,.95,.95);
  legend1->AddEntry(g_dec,"final filter","l");
  legend1->Draw();
  c2->Print(Form("%s/ZAC-dec.pdf",resDir));  
  return 0;
}

