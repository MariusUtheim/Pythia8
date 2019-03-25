// Headers and Namespaces.
#include "Pythia8/Pythia.h" 	// Include Pythia headers.
using namespace Pythia8;	// Let Pythia8:: be implicit.

// Extracting particles from the jet for analysis.
void getJetEvent( Event& event, Event& jetEvent, SlowJet& slowJet, int& i) {
  // Copy over all particles that exist within the jet.
  jetEvent.reset();
  vector<int> constit = slowJet.constituents(i);
  for (int j = 0; j < int(constit.size()); ++j)
    jetEvent.append( event[constit[j]] );
}

// Getting the particles that are not part of the subjet.
// Just a test, not implemented.
/*
  void removeSubjet( Event& jetEvent, SlowJet& slowJet, int& i) {
  vector<int> constit = slowJet.constituents(i);
  for (int j = 0; j < int(constit.size()); ++j)
    jetEvent.remove( constit[j], constit[j] );
}
*/


int main() {
	// Begin main program.
	// Set up generation.
  Pythia pythia;		// Declare Pythia object

  pythia.readString("Beams:eCM = 14000."); // 14 TeV CM energy.
  pythia.readString("PhaseSpace:pTHatMin = 120.");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("WeakBosonAndParton:qqbar2gmZg=on"); // qqbar to gamma Z0 g
  pythia.readString("WeakBosonAndParton:qg2gmZq=on"); // qq to gamma Z0 q
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");
  //pythia.readString("PartonLevel:FSRinResonances = off");
  pythia.init(); 		// Initialize; incoming pp beams is default.

  // Jet level event record.
  Event jetEvent;
  jetEvent.init("Jet level event record", &pythia.particleData);

  // Parameter for the jet finders.
  double etaMax      = 5.;
  double radius      = 1.;
  double pTjetMin    = 200.;
  double subRadius   = 0.3;
  double pTsubjetMin = 10.;
  int nSel           = 2.;
  // Set up SlowJet jet finder, with anti-kT clustering
  // and pion mass assumed for non-photons..
  SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);
  SlowJet slowSubJet( -1, subRadius, pTsubjetMin, etaMax, nSel, 1);


  // Histograms for the Z bosons and the jets they produce.
  Hist pTdiff("pT difference between jets and Z", 100, 0., 100.);
  Hist dist("R distance between jet and Z", 100, 0., 0.1);
  Hist pTsubjets("pT for subjets", 100, 0., 500.);
  Hist pTsubDiff("pT difference between jets and subjets", 100, 0., 400.);
  Hist dist2("R distance between jet and subjet", 100, 0., 1.5);
  Hist num("Number of subjets from jets", 100, 0., 6.);
  Hist dist3("R distance between daughters", 100, 0., 4.);
  Hist dist4("R distance between subjets", 100, 0., 4.);
  Hist dist5("R distance between daughters and subjets", 100, 0., 4.);
  // Number of events to run.
  int runs = 100;
  int idbg = 0;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < runs; ++iEvent){
    if (!pythia.next()) continue;

    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 23) iZ = i;

    double Z_y = pythia.event[iZ].y();
    double Z_phi = pythia.event[iZ].phi();
    double Z_pT = pythia.event[iZ].pT();


    // Analyze SlowJet jet properties.
    slowJet.analyze( pythia.event );

    // Comparing the final Z boson with the jets that have
    // been found and analyzing their properties.
    for (int i = 0; i < slowJet.sizeJet(); ++i) {
      double jet_y = slowJet.y(i);
      double jet_phi = slowJet.phi(i);
      double jet_pT = slowJet.pT(i);
      double dy = jet_y - Z_y;
      double dPhi = abs( jet_phi - Z_phi );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dy) + pow2(dPhi) );

      // If the Z0 and jet are within a certain distance, they are
      // assumed to be correlated and the histograms are filled.
      if (dR < 0.1) {
        pTdiff.fill( Z_pT - jet_pT );
        dist.fill( dR );

        // Finding the daughter particles to the Z.
        int iDau1 = pythia.event[iZ].daughter1();
        int iDau2 = pythia.event[iZ].daughter2();
        double dau1_y = pythia.event[iDau1].y();
        double dau2_y = pythia.event[iDau2].y();
        double dau1_phi = pythia.event[iDau1].phi();
        double dau2_phi = pythia.event[iDau2].phi();
        double dau1_pT = pythia.event[iDau1].pT();
        double dau2_pT = pythia.event[iDau2].pT();

        // Calculating delta R for the daughters.
        double dDau_y = dau1_y - dau2_y;
        double dDau_phi = abs( dau1_phi - dau2_phi );
        if (dDau_phi > M_PI) dDau_phi = 2. * M_PI - dDau_phi;
        double dDau_R = sqrt( pow2(dDau_y) + pow2(dDau_phi) );

        dist3.fill( dDau_R );

        // Analyzing the jet and finding subjets and their properties.
        getJetEvent( pythia.event, jetEvent, slowJet, i );
        slowSubJet.analyze( jetEvent );

        num.fill( slowSubJet.sizeJet() );

        // Fill histograms with subjet information.
        for (int j = 0; j < slowSubJet.sizeJet(); ++j) {
          pTsubjets.fill( slowSubJet.pT(j) );

          // Calculating delta R for jets and subjets.
          double dy2 = jet_y - slowSubJet.y(j);
          double dPhi2 = abs( jet_phi - slowSubJet.phi(j) );
          if (dPhi2 > M_PI) dPhi2 = 2. * M_PI - dPhi2;
          double dR2 = sqrt( pow2(dy2) + pow2(dPhi2) );
          dist2.fill( dR2 );
          pTsubDiff.fill( jet_pT - slowSubJet.pT(j) );

          // Calculating delta R for daughters and subjets.
          double dy4 = dau1_y - slowSubJet.y(j);
          double dPhi4 = abs( dau1_y - slowSubJet.phi(j) );
          if (dPhi4 > M_PI) dPhi4 = 2. * M_PI - dPhi4;
          double dR4 = sqrt( pow2(dy4) + pow2(dPhi4) );
          dist5.fill( dR4 );

          double dy5 = dau2_y - slowSubJet.y(j);
          double dPhi5 = abs( dau2_y - slowSubJet.phi(j) );
          if (dPhi5 > M_PI) dPhi5 = 2. * M_PI - dPhi5;
          double dR5 = sqrt( pow2(dy5) + pow2(dPhi5) );
          dist5.fill( dR5 );

          // Calculating delta R for subjets.
          for (int k = j + 1; k < slowSubJet.sizeJet(); ++k) {
            double dy3 = slowSubJet.y(j) - slowSubJet.y(k);
            double dPhi3 = abs( slowSubJet.phi(j) - slowSubJet.phi(k) );
            if (dPhi3 > M_PI) dPhi3 = 2. * M_PI - dPhi3;
            double dR3 = sqrt( pow2(dy3) + pow2(dPhi3) );
            dist4.fill( dR3 );
          }
        }

        // Debug info.
        if (++idbg < 5) {
          pythia.process.list();
          cout << " Z in " << iZ << " and daughters in " << iDau1
               << " and " << iDau2 << endl;
          pythia.event.list();
          cout << fixed << "   Z_pT = " << Z_pT << "   Z_y = " << Z_y
               << "   Z_phi = " << Z_phi << endl;
          cout << fixed << " jet_pT = " << jet_pT  << " jet_y = " << jet_y
               << " jet_phi = " << jet_phi << endl;
          slowJet.list();
          cout << fixed << " dau1_pT = " << dau1_pT  << " dau1_y = " << dau1_y
               << " dau1_phi = " << dau1_phi << endl;
          cout << fixed << " dau2_pT = " << dau2_pT  << " dau2_y = " << dau2_y
               << " dau2_phi = " << dau2_phi << endl;
          slowSubJet.list();
          for (int j = 0; j < slowSubJet.sizeJet(); ++j) {
            double dy4 = dau1_y - slowSubJet.y(j);
            //double dPhi4 = abs( dau1_y - slowSubJet.phi(j) );
            double dPhi4 = abs( dau1_phi - slowSubJet.phi(j) );
            if (dPhi4 > M_PI) dPhi4 = 2. * M_PI - dPhi4;
            double dR4 = sqrt( pow2(dy4) + pow2(dPhi4) );
            double dy5 = dau2_y - slowSubJet.y(j);
            //double dPhi5 = abs( dau2_y - slowSubJet.phi(j) );
            double dPhi5 = abs( dau2_phi - slowSubJet.phi(j) );
            if (dPhi5 > M_PI) dPhi5 = 2. * M_PI - dPhi5;
            double dR5 = sqrt( pow2(dy5) + pow2(dPhi5) );
            cout << fixed << " jet " << j << " R_1 = " << dR4
                 << " R_2 = " << dR5 << endl;
          }
        }

      }
    }

  // End of event loop. Statistics. Histograms.
  }
  pythia.stat();

  cout << pTdiff << pTsubjets << pTsubDiff << dist << dist2
  << dist3 << dist4 << dist5 << num;

  return 0;
}  // End main program with error-free return.
