// test215.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Space-time evolution of the hadronization process.
// Output plotted by Python/Matplotlib/pyplot.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events.
  int nEvent  = 1000;
  int nAbort  = 5;
  int nStatus = 20;

  // Time/distance steps for histogramming and other counters.
  double tScale[90], sScale[100];
  int    ntPri[90], ntSec[90], ntAll[90], nrPri[90], nrSec[90], nrAll[90],
         nsPri[100], nsSec[100], nsAll[100];
  for (int j = 0; j < 90; ++j) {
    tScale[j] = pow( 10., (j + 0.5) / 6. - 15.);
    ntPri[j] = ntSec[j] = ntAll[j] = nrPri[j] = nrSec[j] = nrAll[j] = 0;
  }
  for (int j = 0; j < 100; ++j) {
    sScale[j] = 0.2 * (j + 0.5);
    nsPri[j] = nsSec[j] = nsAll[j] = 0;
  }
  bool   isPri, isSec;
  int    iDau;
  double tProd, tDec, rProd, rDec, sProd, sDec, zProd, yProd;

  // Histogram particle time and space distributions.
  Hist tPri("primary hadrons",   90, 1e-15, 1., true);
  Hist tSec("secondary hadrons", 90, 1e-15, 1., true);
  Hist tAll("all hadrons",       90, 1e-15, 1., true);
  Hist rPri("primary hadrons",   90, 1e-15, 1., true);
  Hist rSec("secondary hadrons", 90, 1e-15, 1., true);
  Hist rAll("all hadrons",       90, 1e-15, 1., true);
  Hist sPri("primary hadrons",   100, 0., 20.);
  Hist sSec("secondary hadrons", 100, 0., 20.);
  Hist sAll("all hadrons",       100, 0., 20.);
  Hist qPri("primary hadrons",   100, 0., 10.);
  Hist qSec("secondary hadrons", 100, 0., 10.);
  Hist qAll("all hadrons",       100, 0., 10.);
  Hist yPri("primary hadrons",   100, -10., 10.);
  Hist ySec("secondary hadrons", 100, -10., 10.);
  Hist yAll("all hadrons",       100, -10., 10.);

  // Pythia generator. Event record shorthand.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("Beams:eCM = 13000.");

  // Process choice.
  pythia.readString("SoftQCD:inelastic = on");
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("PhaseSpace:pTHatMin = 100.");

  // Vertex selection choices.
  pythia.readString("Fragmentation:setVertices = on");

  // Initialization.
  pythia.init();

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Loop over all particles and analyze the hadrons.
    for (int i = 0; i < pythia.event.size(); ++i)
    if (event[i].isHadron() && event[i].statusAbs() > 20) {
      isPri = (event[i].statusAbs() / 10 == 8);
      isSec = (event[i].statusAbs() / 10 == 9);
      iDau = event[i].daughter1();
      if (!isPri && !isSec && ++nStatus < 20) cout << " status "
        << event[i].status() << " for hadron " << event[i].id()
        << " with mother id " << event[event[i].mother1()].id() << endl;

      // Particle production and decay time and radius, in units of m.
      tProd = 1e-3 * event[i].tProd();
      tDec  = (event[i].status() > 0.) ? 1. : 1e-3 * event[iDau].tProd();
      rProd = 1e-3 * sqrt( pow2(event[i].xProd()) + pow2(event[i].yProd()) );
      rDec  = (event[i].status() > 0.) ? 1.
        : 1e-3 * sqrt( pow2(event[iDau].xProd()) + pow2(event[iDau].yProd()) );
      sProd = 1e15 * rProd;
      sDec  = 1e15 * rDec;
      zProd = 1e-3 * event[i].zProd();
      yProd = 0.5 * log( max(1e-30, tProd + zProd)
            / max(1e-30, tProd - zProd) );

      // Check whether it exists at any of the time slots and count up.
      for (int j = 0; j < 90; ++j) if (tProd < tScale[j] && tDec > tScale[j]) {
        if (isPri) ++ntPri[j];
        if (isSec) ++ntSec[j];
        ++ntAll[j];
      }

      // Check whether it exists at any of the radii slots and count up.
      for (int j = 0; j < 90; ++j) if (rProd < tScale[j] && rDec > tScale[j]) {
        if (isPri) ++nrPri[j];
        if (isSec) ++nrSec[j];
        ++nrAll[j];
      }
      for (int j = 0; j < 100; ++j) if (sProd < sScale[j] && sDec > sScale[j]) {
        if (isPri) ++nsPri[j];
        if (isSec) ++nsSec[j];
        ++nsAll[j];
      }

      // Histogram production r and y_tau.
      if (isPri) qPri.fill( sProd);
      if (isSec) qSec.fill( sProd);
      qAll.fill( sProd);
      if (isPri) yPri.fill( yProd);
      if (isSec) ySec.fill( yProd);
      yAll.fill( yProd);
    }

  // End of event loop. Final statistics.
  }
  pythia.stat();

  // Histogram particle time and space distributions.
  double wt = 1. / nEvent;
  for (int j = 0; j < 90; ++j) {
    tPri.fill( tScale[j], wt * ntPri[j]);
    tSec.fill( tScale[j], wt * ntSec[j]);
    tAll.fill( tScale[j], wt * ntAll[j]);
    rPri.fill( tScale[j], wt * nrPri[j]);
    rSec.fill( tScale[j], wt * nrSec[j]);
    rAll.fill( tScale[j], wt * nrAll[j]);
  }
  for (int j = 0; j < 100; ++j) {
    sPri.fill( sScale[j], wt * nsPri[j]);
    sSec.fill( sScale[j], wt * nsSec[j]);
    sAll.fill( sScale[j], wt * nsAll[j]);
  }
  qPri *= 10. * wt;
  qSec *= 10. * wt;
  qAll *= 10. * wt;
  yPri *= 5. * wt;
  ySec *= 5. * wt;
  yAll *= 5. * wt;
  cout << tPri << tSec << tAll << rPri << rSec << rAll
       << sPri << sSec << sAll << qPri << qSec << qAll
       << yPri << ySec << yAll;

  // Write Python code that can generate a PDF file with the distributions.
  HistPlot hpl("test215plot");
  hpl.frame( "out215plot", "Hadron multiplicity at different times",
    "$ct$ (m)", "$N_{\\mathrm{hadron}}(t)$");
  hpl.add( tPri, "-");
  hpl.add( tSec, "-");
  hpl.add( tAll, "-");
  hpl.plot();
  hpl.frame( "", "Hadron multiplicity at different radii",
    "$r$ (m)", "$N_{\\mathrm{hadron}}(r)$");
  hpl.add( rPri, "-");
  hpl.add( rSec, "-");
  hpl.add( rAll, "-");
  hpl.plot();
  hpl.frame( "", "Hadron multiplicity at different radii",
    "$r$ (fm)", "$N_{\\mathrm{hadron}}(r)$");
  hpl.add( sPri, "-");
  hpl.add( sSec, "-");
  hpl.add( sAll, "-");
  hpl.plot();
  hpl.frame( "", "Hadron production at different radii",
    "$r$ (fm)", "$\\mathrm{d}N_{\\mathrm{hadron}} / \\mathrm{d}r$");
  hpl.add( qPri, "-");
  hpl.add( qSec, "-");
  hpl.add( qAll, "-");
  hpl.plot();
  hpl.frame( "", "Hadron production at different $y_{\\tau}$",
    "$r$ (fm)", "$\\mathrm{d}N_{\\mathrm{hadron}} / \\mathrm{d}y_{\\tau}$");
  hpl.add( yPri, "-");
  hpl.add( ySec, "-");
  hpl.add( yAll, "-");
  hpl.plot();

  // Done.
  return 0;
}
