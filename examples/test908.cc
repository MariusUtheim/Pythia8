// test908.cc
// Jet finding for QCD jets and for top.
// Exercise 8.2.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events. QCD or top.
  int nEvent = 10000;
  bool selectQCD = true;

  // Set up and initialize Pythia.
  Pythia pythia;
  pythia.readString("Beams:Ecm = 13000.");
  if (selectQCD) {
    pythia.readString("HardQCD:all = on");
    pythia.readString("PhaseSpace:pTHatMin = 250.");
  } else {
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
  }
  pythia.readString("Next:numberShowEvent = 0");
  pythia.init();
  Event& event = pythia.event;

  // Jet finders.
  double radius = 0.5;
  double pTJmin = 20.;
  double etaMax = 10.;
  SlowJet akJet( -1, radius, pTJmin, etaMax);
  SlowJet caJet(  0, radius, pTJmin, etaMax);
  SlowJet ktJet(  1, radius, pTJmin, etaMax);

  // Histograms for jets.
  Hist nJetak( "anti-kT",  20, -0.5, 19.5);
  Hist nJetca( "Cam/Aach", 20, -0.5, 19.5);
  Hist nJetkt( "kT",       20, -0.5, 19.5);
  Hist ptJetak( "anti-kT",  100, 0., 1000.);
  Hist ptJetca( "Cam/Aach", 100, 0., 1000.);
  Hist ptJetkt( "kT",       100, 0., 1000.);
  Hist dptJetca( "CA - anti", 100, -50., 50.);
  Hist dptJetka( "kt - anti", 100, -50., 50.);
  Hist dptJetkc( "kt - CA",   100, -50., 50.);

  // Generate events.
  for ( int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if ( !pythia.next() ) continue;

    // Analyze jets in current event: number and hardest.
    akJet.analyze( event);
    caJet.analyze( event);
    ktJet.analyze( event);
    int szak =  akJet.sizeJet();
    int szca =  caJet.sizeJet();
    int szkt =  ktJet.sizeJet();
    double pTak = (szak > 0) ? akJet.pT(0) : 0.;
    double pTca = (szca > 0) ? caJet.pT(0) : 0.;
    double pTkt = (szkt > 0) ? ktJet.pT(0) : 0.;

    // Fill histograms.
    nJetak.fill( szak );
    nJetca.fill( szca );
    nJetkt.fill( szkt );
    ptJetak.fill( pTak );
    ptJetca.fill( pTca );
    ptJetkt.fill( pTkt );
    dptJetca.fill( pTca - pTak );
    dptJetka.fill( pTkt - pTak );
    dptJetkc.fill( pTkt - pTca );

  // End of event loop. Print histograms.
  }
  pythia.stat();
  cout << nJetak << nJetca << nJetkt << ptJetak << ptJetca << ptJetkt
       << dptJetca << dptJetka << dptJetkc;

  // Python code for plotting distributions.
  HistPlot hpl("test908plot");
  hpl.frame( "out908plot", "Number of jets", "$n_{\\mathrm{jet}}$",
    "Probability");
  hpl.add( nJetak );
  hpl.add( nJetca );
  hpl.add( nJetkt );
  hpl.plot();
  hpl.frame( "", "Transverse momentum of the hardest jet",
    "$p_{\\perp\\mathrm{jet}}$ (GeV)", "Probability");
  hpl.add( ptJetak );
  hpl.add( ptJetca );
  hpl.add( ptJetkt );
  hpl.plot();
  hpl.frame( "", "Transverse momentum difference of the hardest jet",
    "$\\Delta p_{\\perp\\mathrm{jet}}$ (GeV)", "Probability");
  hpl.add( dptJetca );
  hpl.add( dptJetka );
  hpl.add( dptJetkc );
  hpl.plot(true);

  // Done.
  return 0;
}
