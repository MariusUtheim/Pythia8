// test225.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// The fusion of two jets into one.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Simple method to do the filling of partons into the event record.

void fillPartons(int type, double ee, double theta, Event& event) {

  // Reset event record to allow for new event.
  event.reset();

  // Information on a g g system, to be hadronized.
  if (type == 1) {
    event.append( 21, 23, 101, 102, 0., 0.,  2. * ee, 2. * ee);
    event.append( 21, 23, 102, 101, 0., 0., -2. * ee, 2. * ee);

    // Information on two g g system, to be hadronized.
  } else if (type == 2) {
    event.append( 21, 23, 101, 102,  ee * sin(theta), 0.,  ee * cos(theta), ee);
    event.append( 21, 23, 102, 101, 0., 0., -ee, ee);
    event.append( 21, 23, 103, 104, -ee * sin(theta), 0.,  ee * cos(theta), ee);
    event.append( 21, 23, 104, 103, 0., 0., -ee, ee);

  // Information on a g g g system, to be hadronized.
  } else if (type == 3) {
    event.append( 21, 23, 101, 102,  ee * sin(theta), 0.,  ee * cos(theta), ee);
    event.append( 21, 23, 102, 103, -ee * sin(theta), 0.,  ee * cos(theta), ee);
    event.append( 21, 23, 103, 101, 0., 0., -2. * ee, 2. * ee);
  }
}

//==========================================================================

int main() {

  // Number of events to generate. Conversion to radians.
  int nEvent = 1000000;
  double eJet = 10.;

  // Generator; shorthand for event and particleData.
  Pythia pythia;
  Event& event      = pythia.event;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize.
  pythia.init();

  // Book histograms. Loop over five event cases.
  Hist Hfx[5], Htheta[5];
  for (int mode = 0; mode < 5; ++mode) {
    Hfx[mode].book( " ", 50, 0., 1.);
    Htheta[mode].book( "  ", 50, 0., 0.5 * M_PI);

    // Begin of event loop.
    int nAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Set up parton system according to case.
      if       (mode == 0) fillPartons( 2, eJet, 0.4, event);
      else if  (mode == 1) fillPartons( 3, eJet, 0.4, event);
      else if  (mode == 2) fillPartons( 3, eJet, 0.2, event);
      else if  (mode == 3) fillPartons( 3, eJet, 0.1, event);
      else if  (mode == 4) fillPartons( 1, eJet, 0.0, event);

      // Generate events. Quit if failure.
      if (!pythia.next() && ++nAbort > 50) {
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }

      // List first few events.
      if (iEvent < 1) event.list();

      // Loop over all final particles in forward hemisphere.
      for (int i = 0; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].pz() > 0.) {

        // Longitudinal momentum and particle flow in event plane.
        double xx = 0.5 * event[i].pz() / eJet;
        double thetaXZ = abs(event[i].thetaXZ());
        if (event[i].isCharged()) Hfx[mode].fill( xx );
        Htheta[mode].fill(thetaXZ, event[i].e() );
      }

    // End of event loop. Normalize histograms. End loop over cases.
    }
    pythia.stat();
    Hfx[mode] *= 50. / nEvent;
    Htheta[mode] *= 100. / (M_PI * nEvent);
  }

  // Plot histograms.
  HistPlot hpl("test225plot");
  hpl.frame( "out225plot", "Charged particle fragmentation function",
    "$x$ (energy fraction)", "$\\mathrm{d}N_{\\mathrm{ch}} / \\mathrm{d}x$");
  hpl.add( Hfx[0], "--,black", "two unconnected g, $E = 10$, $\\theta = 0.4$");
  hpl.add( Hfx[1], "-", "two connected g, $E = 10$, $\\theta = 0.4$");
  hpl.add( Hfx[2], "-", "two connected g, $E = 10$, $\\theta = 0.2$");
  hpl.add( Hfx[3], "-", "two connected g, $E = 10$, $\\theta = 0.1$");
  hpl.add( Hfx[4], "-,black", "one g, $E = 20$");
  hpl.plot(true);
  hpl.frame( "", "Energy flow in half event plane",
    "$\\theta$", "$\\mathrm{d}E / \\mathrm{d}\\theta$");
  hpl.add( Htheta[0], "--,black",
    "two unconnected g, $E = 10$, $\\theta = 0.4$");
  hpl.add( Htheta[1], "-",
    "two connected g, $E = 10$, $\\theta = 0.4$");
  hpl.add( Htheta[2], "-",
    "two connected g, $E = 10$, $\\theta = 0.2$");
  hpl.add( Htheta[3], "-",
    "two connected g, $E = 10$, $\\theta = 0.1$");
  hpl.add( Htheta[4], "-,black", "one g, $E = 20$");
  hpl.plot(true);

  // Done.
  return 0;
}
