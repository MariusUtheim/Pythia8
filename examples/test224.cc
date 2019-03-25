// test224.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Study CR effects.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent  = 1000000;
  int nAbort  = 5;

  // Book histograms.
  Hist HnMPI[2], Hnch[2], HpT[2], Hjet[2], Hnchr[2];

  // Loop over runs with or without MPI. Book histograms.
  for (int mode = 0; mode < 2; ++mode) {
    HnMPI[mode].book( "  ", 24, 0.5, 24.5);
    Hnch[mode].book( "  ", 24, 0.5, 24.5);
    Hnchr[mode].book( "  ", 99, 0.5, 792.5);
    HpT[mode].book( "  ", 50, 0., 5.);
    Hjet[mode].book( "  ", 24, 0.5, 24.5);

    // Generator setup.
    Pythia pythia;
    Event& event = pythia.event;
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("SoftQCD:nonDiffractive = on");
    if (mode == 0) pythia.readString("ColourReconnection:reconnect = off");
    pythia.init();

    // Anti-kT jet finder.
    double etaMax   = 5.;
    double radius   = 0.6;
    double pTjetMin = 20.;
    int    nSel     = 2;
    SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

    // Begin event loop.
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Quit if too many failures.
      if (!pythia.next()) {
        if (++iAbort < nAbort) continue;
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }

      // Number of MPI's.
      int nMPI = pythia.info.nMPI();
      HnMPI[mode].fill( nMPI );

      // Charged multiplicity and charged pT spectrum.
      int nch = 0;
      for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) {
        ++nch;
        HpT[mode].fill( event[i].pT() );
      }
      Hnch[mode].fill( nMPI, nch );
      Hnchr[mode].fill( nch );

      // Number of jets.
      slowJet. analyze( pythia.event );
      Hjet[mode].fill( nMPI, slowJet.sizeJet() );
    }
    // End of event loop. Normalize histograms.
    pythia.stat();
    Hnch[mode] /= HnMPI[mode];
    Hjet[mode] /= HnMPI[mode];
    HnMPI[mode] /= nEvent;
    Hnchr[mode] *= 0.25 / nEvent;
    HpT[mode] *= 10. / nEvent;
  }

  // Plot histograms.
  HistPlot hpl("test224plot");
  hpl.frame( "out224plot", "Number of MPIs at 13 TeV",
    "$n_{\\mathrm{MPI}}$", "$\\mathrm{d}P / \\mathrm{d}n_{\\mathrm{MPI}}$");
  hpl.add( HnMPI[0], "-", "CR off");
  hpl.add( HnMPI[1], "-,red", "CR on");
  hpl.plot(true);
  hpl.frame( "", "Charged multiplicity distribution at 13 TeV",
    "$n_{\\mathrm{ch}}$", "$\\mathrm{d}P / \\mathrm{d}n_{\\mathrm{ch}}$");
  hpl.add( Hnchr[0], "-", "CR off");
  hpl.add( Hnchr[1], "-,red", "CR on");
  hpl.plot(true);
  hpl.frame( "", "Charged multiplicity dependence on number of MPIs at 13 TeV",
    "$n_{\\mathrm{MPI}}$", "$\\langle n_{\\mathrm{ch}} \\rangle$");
  hpl.add( Hnch[0], "-", "CR off");
  hpl.add( Hnch[1], "-,red", "CR on");
  hpl.plot();
  hpl.frame( "", "Charged-particle $p_{\\perp}$ spectrum at 13 TeV",
    "$p_{\\perp}$ (GeV)", "$\\mathrm{d}N / \\mathrm{d} p_{\\perp}$");
  hpl.add( HpT[0], "-", "CR off");
  hpl.add( HpT[1], "-,red", "CR on");
  hpl.plot(true);
  hpl.frame( "", "Jet multiplicity dependence on number of MPIs at 13 TeV",
    "$n_{\\mathrm{MPI}}$", "$\\langle n_{\\mathrm{jet}} \\rangle$"
    " (anti-$k_{\\perp}, R < 0.6, p_{\\perp} > 20$");
  hpl.add( Hjet[0], "-", "CR off");
  hpl.add( Hjet[1], "-,red", "CR on");
  hpl.plot();

  // Done.
  return 0;
}
