// File: hist.cc
// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.
// Copyright (C) 2012 Torbjorn Sjostrand

// Stdlib header file for input and output.
#include <iostream>
#include <iomanip>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
//#include "TH1.h"

// ROOT, for interactive graphics.
//#include "TVirtualPad.h"
//#include "TApplication.h"

// ROOT, for saving file.
//#include "TFile.h"


using namespace Pythia8;

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  //TApplication theApp("hist", &argc, argv);

  Pythia pythia;

  pythia.readString("Beams:eCM = 13000.");

  pythia.readString("HiggsSM:gmgm2H = on");
  pythia.readString("25:onMode = off");
  pythia.readString("25:onIfMatch = 22 22");

  // Force heavy Higgs to be narrow.
  pythia.readString("25:m0 = 750.");
  pythia.readString("25:mWidth = 1.");
  pythia.readString("25:doForceWidth = on");




  // Set up photon beams with pointlike flux.
  PDF* pdfAPtr = new ProtonPoint ( 2212);
  //PDF* pdfAPtr  = new NNPDF ( 2212, 4, "../share/Pythia8/xmldoc/", 0);
  //PDF* pdfAPtr  = new LHAGrid1 ( 2212, "CT14qed_proton_0000.dat", "../share/Pythia8/xmldoc/", 0);

  //PDF* pdfBPtr = new ProtonPoint ( 2212);
  //PDF* pdfBPtr = new NNPDF ( 2212, 4, "../share/Pythia8/xmldoc/", 0);
  PDF* pdfBPtr  = new LHAGrid1 ( 2212, "CT14qed_proton_0000.dat", "../share/Pythia8/xmldoc/", 0);

  pythia.setPDFPtr ( pdfAPtr , pdfBPtr );

  pythia.readString("PartonLevel:MPI = off");
  // Do not allow FSR if you want to avoid gamma -> f fbar
  // in Higgs decay.
  pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("PartonLevel:ISR = off");
  // The kT spectrum you get is not the right one, so misleading.
  //pythia.readString("BeamRemnants:primordialKT = off");


  //pythia.readString("BeamRemnants:primordialKTsoft = 0.1");

  // Switch to force either or both remnants to be original proton.
  pythia.settings.addMode("BeamRemnants:unresolvedHadron",
    0, true, true, 0, 3);
  // Force unresolved: 0 = none, 1 = side 1, 2 = side 2, 3 = both.
  pythia.readString("BeamRemnants:unresolvedHadron = 1");


  // No event printout.
  pythia.readString("Next:numberShowEvent = 5");

  pythia.init();

  //TH1F *pt_Pythia8   		= new TH1F("pt_Pythia8"  ,"", 100 ,0. ,10.);
  Hist pt_Pythia8("pt_Pythia8"  , 100 ,0. ,10.);

  // Begin event loop. Generate event; skip if generation aborted.

  int nevt=5000;
  bool firstEvent = true;
  for (int iEvent = 0; iEvent < nevt; ++iEvent) {

    if (!pythia.next()) continue;

    float pt_tot = 0;

    for (int i = 0; i < pythia.event.size(); ++i){


	if(pythia.event[i].id() == 25 && pythia.event[i].status()==-62){
		pt_tot = pythia.event[i].pT();
	}

    }

    pt_Pythia8.fill( pt_tot );

  }

  // Statistics on event generation.
  pythia.stat();

  cout << pt_Pythia8;

  //std::cout << "\nDouble click on the histogram window to quit.\n";
  //gPad->WaitPrimitive();

  // Done.
  return 0;
}
