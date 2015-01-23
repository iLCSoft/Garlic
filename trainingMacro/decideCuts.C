#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TClass.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDatime.h"
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

using std::cout;
using std::endl;

void decideCuts() {

  TCanvas* cc = new TCanvas();
  TString plname = "decideCuts.ps";
  cc->Print(plname+"[");


  const int MAXPDG=2;
  const int MAXV=20;
  const int MAXE=30;
  const int MAXR=4;

  float evtFracs[MAXR] = {0.01, 0.05, 0.95, 0.99};

  std::map < std::pair < int, std::string > , std::vector < std::pair < float, std::vector < float > > > > allInfo;

  for (int jpdg=0; jpdg<2; jpdg++) {

    int ipdg(0);
    std::vector < std::string > gamma_files;

    if ( jpdg==0 ) {
      ipdg=22;
      gamma_files.push_back("clusterParams_gamma_0.2.root");
      gamma_files.push_back("clusterParams_gamma_0.3.root");
      gamma_files.push_back("clusterParams_gamma_0.4.root");
      gamma_files.push_back("clusterParams_gamma_0.5.root");
      gamma_files.push_back("clusterParams_gamma_1.root");
      gamma_files.push_back("clusterParams_gamma_2.root");
      gamma_files.push_back("clusterParams_gamma_3.root");
      gamma_files.push_back("clusterParams_gamma_4.root");
      gamma_files.push_back("clusterParams_gamma_5.root");
      gamma_files.push_back("clusterParams_gamma_10.root");
      gamma_files.push_back("clusterParams_gamma_15.root");
      gamma_files.push_back("clusterParams_gamma_20.root");
      gamma_files.push_back("clusterParams_gamma_25.root");
      gamma_files.push_back("clusterParams_gamma_30.root");
      gamma_files.push_back("clusterParams_gamma_35.root");
      gamma_files.push_back("clusterParams_gamma_40.root");
      gamma_files.push_back("clusterParams_gamma_50.root");
      gamma_files.push_back("clusterParams_gamma_60.root");
      gamma_files.push_back("clusterParams_gamma_70.root");
      gamma_files.push_back("clusterParams_gamma_80.root");
      gamma_files.push_back("clusterParams_gamma_90.root");
      gamma_files.push_back("clusterParams_gamma_100.root");
    } else if ( jpdg==1 ) {
      ipdg=11;
      gamma_files.push_back("clusterParams_e-_0.2.root");
      gamma_files.push_back("clusterParams_e-_0.3.root");
      gamma_files.push_back("clusterParams_e-_0.4.root");
      gamma_files.push_back("clusterParams_e-_0.5.root");
      gamma_files.push_back("clusterParams_e-_1.root");
      gamma_files.push_back("clusterParams_e-_2.root");
      gamma_files.push_back("clusterParams_e-_3.root");
      gamma_files.push_back("clusterParams_e-_4.root");
      gamma_files.push_back("clusterParams_e-_5.root");
      gamma_files.push_back("clusterParams_e-_10.root");
      gamma_files.push_back("clusterParams_e-_15.root");
      gamma_files.push_back("clusterParams_e-_20.root");
      gamma_files.push_back("clusterParams_e-_25.root");
      gamma_files.push_back("clusterParams_e-_30.root");
      gamma_files.push_back("clusterParams_e-_35.root");
      gamma_files.push_back("clusterParams_e-_40.root");
      gamma_files.push_back("clusterParams_e-_50.root");
      gamma_files.push_back("clusterParams_e-_60.root");
      gamma_files.push_back("clusterParams_e-_70.root");
      gamma_files.push_back("clusterParams_e-_80.root");
      gamma_files.push_back("clusterParams_e-_90.root");
      gamma_files.push_back("clusterParams_e-_100.root");

    }

    if ( int(gamma_files.size()) > MAXE ) {
      cout << "increase MAXE to at least " << gamma_files.size() << endl;
      return;
    }

    int   npt[MAXPDG]={0};


    for (size_t ifile=0; ifile<gamma_files.size(); ifile++) {
      npt[jpdg]++;
      cout << "file name = " << gamma_files[ifile] << endl;
      TFile* fin = new TFile( gamma_files[ifile].c_str(),"read" );
      int minbin(-1);
      int maxbin(-1);
      int ivar(0);
      float thisEnergy(0);
      bool firstOne=true;
      TIter next(fin->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH2F")) continue;
        TH2F *h = dynamic_cast <TH2F*> (key->ReadObj());
        //cout << h->GetName() << " " << h->GetEntries() << endl;
        TString nn = h->GetName();
        TObjArray* toa = nn.Tokenize("_");
        TString vname = ((TObjString*) toa->At(1))->GetString();
        if ( toa->GetEntries() > 2 ) {
          for (int i=2; i<toa->GetEntries(); i++) {
            vname+="_";
            vname+=((TObjString*) toa->At(i))->GetString();
          }
        }
        std::string al( vname );
        std::pair < int, std::string > thisOne( ipdg, al );
        if ( allInfo.find(thisOne)==allInfo.end() ) {
          std::vector < std::pair < float, std::vector < float > > > dd;
          allInfo[thisOne] = dd;
        }
        if ( firstOne ) {
          // decide energy cuts
          TH1D* hd = h->ProjectionX();
          float tot = hd->Integral(1, hd->GetNbinsX());
          float integ(0);
          for (int ib=1; ib<hd->GetNbinsX(); ib++) {
            integ+=hd->GetBinContent(ib)/tot;
            if ( minbin<0 && integ>0.05 ) minbin=ib;
            if ( maxbin<0 && integ>0.95 ) maxbin=ib;
          }
          //cout << "min, max bins, values = " << minbin << " " << maxbin << " , " << hd->GetBinCenter(minbin) << " " << hd->GetBinCenter(maxbin) << endl;

          //      hd->GetXaxis()->SetRangeUser( hd->GetBinCenter(minbin), hd->GetBinCenter(maxbin) );
	  //          energies[ipdg][ifile]=hd->GetMean();

	  thisEnergy=hd->GetMean();

          //cout << "mean energy = " << thisEnergy << endl;

          firstOne=false;
          delete hd;
        } // first one

        TH1D* hy = h->ProjectionY("_py", minbin, maxbin);

        float tot = hy->Integral(1, hy->GetNbinsX() );

        float evtFracPoints[MAXR]; for (int i=0; i<MAXR; i++) evtFracPoints[i]=999;

        float lower(0);
        for (int i=1; i<=hy->GetNbinsX(); i++) {
          lower+=hy->GetBinContent(i)/tot;
          for (int jj=0; jj<MAXR; jj++) {
            if ( evtFracPoints[jj]>998 && lower > evtFracs[jj] ) {
              if ( evtFracs[jj]<0.5 ) {
                evtFracPoints[jj]=hy->GetBinLowEdge(i) - hy->GetBinWidth(i);
              } else {
                evtFracPoints[jj]=hy->GetBinLowEdge(i) + 2*hy->GetBinWidth(i);
              }
            }
          }
        }

        std::vector < float > percs;
        for (int jj=0; jj<MAXR; jj++) percs.push_back( evtFracPoints[jj] );

        allInfo[thisOne].push_back( std::pair < float, std::vector < float > > (thisEnergy,percs) );

        delete h;
        ivar++;
      }
      fin->Close();
    }

  }

  // got all the information, now fit

  TDatime date;

  ofstream myfile;
  myfile.open ("garlicCuts.txt");
  myfile << "#####################################################" << endl;
  myfile << "#### photon and electron cluster cuts for GARLIC ####" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "## these are the ranges within which we expect     ##" << endl;
  myfile << "## various cluster parameters to fall for true     ##" << endl;
  myfile << "## electromagnetic clusters.                       ##" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "## 1,5,95,99 percentiles are given:                ##" << endl;
  myfile << "## the values below which that percentage of       ##" << endl;
  myfile << "## EM clusters lie.                                ##" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "## These percentiles are modelled as polynomial    ##" << endl;
  myfile << "## functions of log10(cluster energy).             ##" << endl;
  myfile << "## The a_i given below are the parameters of these ##" << endl;
  myfile << "## polynomials. In principle, any degree of        ##" << endl;
  myfile << "## poly can be used (in practice, usually pol3     ##" << endl;
  myfile << "## is a good choice)                               ##" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "## f(E)=a0+a1*log10(E)+a2*pow(log10(E),2)+ ....    ##" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "## PDG=11 (electron) 22 (photon)                   ##" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "## PDG var_name percentile a0 a1 a2 a3...          ##" << endl;
  myfile << "##                                                 ##" << endl;
  myfile << "#####################################################" << endl;
  myfile << "## this file created on " << date.AsString() << endl;
  myfile << "#####################################################" << endl;

  int ipp(-1);

  for ( std::map < std::pair < int, std::string > , std::vector < std::pair < float, std::vector < float > > > > ::iterator itt = allInfo.begin(); itt!=allInfo.end(); itt++) {

    if ( itt->first.first != ipp ) {
      ipp = itt->first.first;
      myfile << "PDG " << ipp << endl;
      myfile << " CLASS DUMMY" << endl;
    }

    cout << "PDG_VAR " << itt->first.first << " " << itt->first.second << endl;

    myfile << "  VAR " << itt->first.second << endl;

    std::vector < std::pair < float, std::vector < float > > > aa = itt->second;
    float value[MAXR][MAXE];
    float logen[MAXE];
    if ( int(aa.size()) > MAXE || int(aa[0].second.size()) > MAXR ) cout << "increase MAXR, MAXE" << endl;
    float pmin( 999999);
    float pmax(-999999);
    for ( size_t j=0; j<aa.size(); j++ ) { // energies
      logen[j] = log10(aa[j].first);
      std::vector < float > bb = aa[j].second;
      for (size_t k=0; k<bb.size(); k++) { // percentiles
        pmin=std::min( pmin, bb[k]);
        pmax=std::max( pmax, bb[k]);
        value[k][j] = bb[k];
      }
    }
    float delta = (pmax-pmin)/5;
    TString name; name+=itt->first.first; name+="_"; name+=itt->first.second;
    TH2F* hsum = new TH2F(name, name,10,-1,2,10,pmin-delta,pmax+delta);
    cc->Clear();
    hsum->Draw();
    for (size_t k=0; k<4; k++) { // percentiles
      // myfile << itt->first.first << " " << itt->first.second << " " << evtFracs[k] << " ";
      myfile << "   EVTFRAC_PARS " << evtFracs[k] << " : ";
      TGraph* grr = new TGraph( aa.size(), logen, value[k]);
      grr->Fit("pol3","q");
      TF1* ff = grr->GetFunction("pol3");
      for (int l=0; l<4; l++)
        myfile << ff->GetParameter(l) << " ";
      myfile << endl;
      grr->SetLineColor(k+1);
      grr->SetMarkerColor(k+1);
      ff->SetLineColor(k+1);
      grr->Draw("plsame");
    }
    cc->Print(plname);
  }

  myfile.close();

  cc->Print(plname+"]");


}
