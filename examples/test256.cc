// test256.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Development of the LowEnergyHadHad class for low-energy hadron-hadron
// collisions, as needed e.g. in rescattering descriptions.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events to generate. Max number of errors. 
  int nEvent = 1000;
  int nAbort = 5;
  int nList  = 1;

  // Minimum invariant mass excess for rescattering.
  double mRescMin = 0.5;

  // Set up inelastic nondiffractive LHC collisions.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();

  // Histograms.
  Hist nResc( "Number of rescatterings in an event", 100, -0.5, 199.5); 
  Hist mExc( "Mass excess in scattering", 100, 0., 50.); 
  Hist nHad( "Number of final-state hadrons given mass excess", 100, 0., 50.); 
  Hist mExcS( "Mass excess in successful scattering", 100, 0., 50.); 
  Hist nHadS( "Number of final-state hadrons when successful", 100, 0., 50.); 
  Hist mExc1( "Mass excess in nondiffractive", 100, 0., 50.); 
  Hist nHad1( "Number of hadrons in nondiffractive", 100, 0., 50.); 
  Hist mExc2( "Mass excess in elastic", 100, 0., 50.); 
  Hist nHad2( "Number of hadrons in elastic", 100, 0., 50.); 
  Hist mExc34( "Mass excess in single diffractive", 100, 0., 50.); 
  Hist nHad34( "Number of hadrons in single diffractive", 100, 0., 50.); 
  Hist mExc5( "Mass excess in double diffractive", 100, 0., 50.); 
  Hist nHad5( "Number of hadrons in double diffractive", 100, 0., 50.); 
  Hist mExc6( "Mass excess in annihilation", 100, 0., 50.); 
  Hist nHad6( "Number of hadrons in annihilation", 100, 0., 50.); 

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

   // Find an original hadron pair with invariant mass above threshold.
   int ncollide = 0;
   int sizeOrig = event.size();
   int sizeBef, sizeAdd;
   if (iEvent < nList) cout << " original event size is " << sizeOrig << endl;
   for (int i1 = 0; i1 < sizeOrig - 1; ++i1) 
   if (event[i1].isFinal() && event[i1].isHadron()) {
     for (int i2 = i1 + 1; i2 < sizeOrig; ++i2) 
     if (event[i2].isFinal() && event[i2].isHadron()) {
       double mExcess = (event[i1].p() + event[i2].p()).mCalc() 
         - event[i1].m() - event[i2].m(); 
       if (mExcess > mRescMin) {
  
         // Pick type of collision.
         int type = 1. + 6. * pythia.rndm.flat();
         //int type = 6;
         if (iEvent < nList) cout << "\n collision between " << i1 << " and " 
            << i2 << " is of type " << type << endl;

         // Perform collision.
         sizeBef = event.size();
         pythia.doLowEnergyHadHad( i1, i2, type);
         sizeAdd = event.size() - sizeBef;

         // Statistics. 
         ++ncollide;
         mExc.fill( mExcess);
         nHad.fill( mExcess, sizeAdd);
         if (sizeAdd > 0) { 
           mExcS.fill( mExcess);
           nHadS.fill( mExcess, sizeAdd);
           if (type == 1) {
             mExc1.fill( mExcess);
             nHad1.fill( mExcess, sizeAdd);
           } else if (type == 2) {
             mExc2.fill( mExcess);
             nHad2.fill( mExcess, sizeAdd);
           } else if (type == 3 || type == 4) {
             mExc34.fill( mExcess);
             nHad34.fill( mExcess, sizeAdd);
           } else if (type == 5) {
             mExc5.fill( mExcess);
             nHad5.fill( mExcess, sizeAdd);
           } else if (type == 6) {
             mExc6.fill( mExcess);
             nHad6.fill( mExcess, sizeAdd);
           }              
         } 
         break;
       }
     }
   } 

   // Check and printout.
   Vec4 pSum;
   for (int i = 1; i < event.size(); ++i) if (event[i].isFinal())
     pSum += event[i].p();
   Vec4 pRef = event[0].p();
   double diff = abs(pSum.px() - pRef.px()) + abs(pSum.py() - pRef.py())
               + abs(pSum.pz() - pRef.pz()) + abs(pSum.e() - pRef.e());
   bool hasError = (diff > 0.1);
   if (iEvent < nList || hasError) {
     cout << " event contained " << ncollide << " collisions " << endl;
     if (hasError) cout << " error, pSum = " << pSum;
     event.list();
   }
   nResc.fill( ncollide);

  // End event loop. Final statistics.
  }
  pythia.stat();
  nHad /= mExc;
  nHadS /= mExcS;
  nHad1 /= mExc1;
  nHad2 /= mExc2;
  nHad34 /= mExc34;
  nHad5 /= mExc5;
  nHad6 /= mExc6;
  cout << nResc << mExc << nHad << mExcS << nHadS << mExc1 << nHad1 
       << mExc2 << nHad2 << mExc34 << nHad34 << mExc5 << nHad5 
       << mExc6 << nHad6;

  // Done.
  return 0;
}
