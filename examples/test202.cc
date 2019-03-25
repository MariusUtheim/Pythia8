// test202.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test program of diffractive scattering in various frameworks.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

// Constants.
const double MPROTON    = 0.9382720;

// Global histograms.
Hist nTries( "number of tries", 100, -0.5, 999.5);
Hist lgxiAll( "log10(xi) all tries", 100, -10., 0.);
Hist lgxiAcc( "log10(xi) accepted fraction", 100, -10., 0.);
Hist tAbsAll( "|t| all tries", 100, 0., 5.);
Hist tAbsAcc( "|t| accepted fraction", 100, 0., 5.);

//==========================================================================

int main() {

  // Choose pp/ppbar, CM energy and number of events.
  bool   ispbarp     = false;
  double eCM         = 10000.;
  int    nEvent      = 10000;

  // Choose processes to generate.
  bool   hasEl       = false;
  bool   hasSD       = false;
  bool   hasDD       = false;
  bool   hasCD       = true;

  // Switch for all elastic scattering scenarios.
  bool   hasCoulomb  = false;

  // TotEl & Diff: 0 = own; 1 = SaS/DL, 2 = MBR; 3 = ABMST; 4 = RPP; -1 = all.
  int    iModeTotEl  = 3;
  int    iModeDiff   = -1;

  // Choices for the "set own" framework.
  bool   OwndampenGap = true;

  // Choices and parameters for the ABMST-based framework.
  // ABMSTmodeSD: 0 = original cut at |t| < 4; 1 = modified with exponential.
  int    ABMSTmodeSD = 1;
  // ABMSTmodeDD: 0 = SD * SD / El; 1 += elastic no dip;
  // 2 += rescale mult * (s/m_p^2)^pow.
  int    ABMSTmodeDD = 1;
  double ABMSTmultDD = 1.;
  double ABMSTpowDD  = 0.1;
  // ABMSTmodeCD: 0 = SD * SD / Tot; 1 += rescale mult * (s/m_p^2)^pow.
  int    ABMSTmodeCD = 0;
  double ABMSTmultCD = 1.;
  double ABMSTpowCD  = 0.1;
  // Damping of small rapidity gaps in all ABMST diffractive processes.
  bool   ABMSTdampenGap = true;
  double ABMSTygap   = 2.;
  double ABMSTypow   = 5.;
  // Minimal slopes for ABMST diffractive scattering.
  bool   ABMSTuseBMin = true;
  double ABMSTbMinSD = 2.;
  double ABMSTbMinDD = 2.;
  double ABMSTbMinCD = 2.;

  // Book histograms.
  Hist tLog("log10(|t|) spectrum",               80,  -6., 2.);
  Hist xiMaxLog( "log10(xi_Max) spectrum",      100, -10., 0.);
  Hist xiMinLog( "log10(xi_Min) spectrum",      100, -10., 0.);
  Hist xiProdLog("log10(xi_1 * xi_2) spectrum", 100, -10., 0.);
  Hist yGap("largest rapidity gap size",        100, -5., 20.);
  Hist yGapMin("smallest rapidity gap size",    100, -5., 20.);

  // Loop over different fluxes for user-set diffractive scenario.
  // Also other options accessible.
  int nFlux = (iModeDiff == 0) ? 7 : 1;
  if (iModeDiff == -1) nFlux = 10;
  for (int iFlux = 1; iFlux <= nFlux; ++iFlux) {
    if (nFlux > 1) cout << "\n Begin iFlux option " << iFlux << endl;
    int iModeDiffNow = iModeDiff;
    if (iModeDiff == -1) iModeDiffNow = (iFlux <= 7) ? 0 : iFlux - 7;
    int iFluxNow = (iModeDiff > 0 || iFlux > 7) ? 1 : iFlux;

    // Generator without header printout. Hard process record shorthand.
    Pythia pythia("../share/Pythia8/xmldoc", false);
    Event& proc = pythia.process;

    // Set up run: total and elastic cross section.
    if (hasEl) pythia.readString("SoftQCD:elastic = on");
    pythia.settings.mode("SigmaTotal:mode", iModeTotEl);
    pythia.settings.flag("SigmaElastic:Coulomb", hasCoulomb);
    if (hasCoulomb) pythia.readString("SigmaElastic:tAbsMin = 4e-5");

    // Set up run: diffractive cross sections.
    if (hasSD) pythia.readString("SoftQCD:singleDiffractive = on");
    if (hasDD) pythia.readString("SoftQCD:doubleDiffractive = on");
    if (hasCD) pythia.readString("SoftQCD:centralDiffractive = on");
    if (hasCD) pythia.readString("SigmaTotal:zeroAXB = off");
    pythia.settings.mode("SigmaDiffractive:mode", iModeDiffNow);
    pythia.settings.mode("SigmaDiffractive:PomFlux", iFluxNow);

    // Set up run: settings for the "set own" scenario.
    pythia.settings.flag("SigmaDiffractive:OwndampenGap", OwndampenGap);

    // Set up run: settings for the ABMST-based scenario.
    pythia.settings.mode("SigmaDiffractive:ABMSTmodeSD", ABMSTmodeSD);
    pythia.settings.mode("SigmaDiffractive:ABMSTmodeDD", ABMSTmodeDD);
    pythia.settings.mode("SigmaDiffractive:ABMSTmodeCD", ABMSTmodeCD);
    pythia.settings.parm("SigmaDiffractive:ABMSTmultDD", ABMSTmultDD);
    pythia.settings.parm("SigmaDiffractive:ABMSTpowDD",  ABMSTpowDD);
    pythia.settings.parm("SigmaDiffractive:ABMSTmultCD", ABMSTmultCD);
    pythia.settings.parm("SigmaDiffractive:ABMSTpowCD",  ABMSTpowCD);
    pythia.settings.flag("SigmaDiffractive:ABMSTdampenGap", ABMSTdampenGap);
    pythia.settings.parm("SigmaDiffractive:ABMSTygap",   ABMSTygap);
    pythia.settings.parm("SigmaDiffractive:ABMSTypow",   ABMSTypow);
    pythia.settings.flag("SigmaDiffractive:ABMSTuseBMin", ABMSTuseBMin);
    pythia.settings.parm("SigmaDiffractive:ABMSTbMinSD", ABMSTbMinSD);
    pythia.settings.parm("SigmaDiffractive:ABMSTbMinDD", ABMSTbMinDD);
    pythia.settings.parm("SigmaDiffractive:ABMSTbMinCD", ABMSTbMinCD);

    // Final common setup. Initialization.
    if (ispbarp) pythia.readString("Beams:idA = -2212");
    pythia.settings.parm("Beams:eCM", eCM);
    pythia.readString("Next:numberCount = 1000000");
    pythia.readString("PartonLevel:all = off");
    int nAbort = 5;
    pythia.init();

    // Begin event loop.
    int iAbort = 0;
    int nxi005 = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events. Quit if too many failures.
      if (!pythia.next()) {
        if (++iAbort < nAbort) continue;
        cout << " Event generation aborted prematurely, owing to error!\n";
        break;
      }

      // Study t distributions of elastic/diffractive events.
      int code = pythia.info.code();
      double tAbsL = log10(abs(pythia.info.tHat()));
      tLog.fill( tAbsL);
      if (code == 106) {
        double uAbsL = log10(abs(pythia.info.uHat()));
        tLog.fill( uAbsL);
      }

      // Study log10(xi) and rapidity gap distributions: SD, DD, CD.
      double xi, xi1, xi2, expGap, expGap3, expGap4, y5;
      if (code == 103 || code == 104) {
        xi  = pow2( max( proc[3].m(), proc[4].m()) / eCM );
        xiMaxLog.fill( log10(xi) );
        yGap.fill( -log(xi) );
        if (xi < 0.05) ++nxi005;
      } else if (code == 105) {
        xi1  = pow2( proc[3].m() / eCM );
        xi2  = pow2( proc[4].m() / eCM );
        xiMaxLog.fill( log10( max(xi1, xi2) ) );
        xiMinLog.fill( log10( min(xi1, xi2) ) );
        expGap = pow2( eCM * MPROTON / (proc[3].m() * proc[4].m()) );
        yGap.fill( log(expGap) );
      } else if (code == 106) {
        xi = pow2( proc[5].m() / eCM);
        xiProdLog.fill( log10(xi) );
        y5   = proc[5].y();
        xi1 = sqrt(xi) * exp( abs(y5));
        xi2 = sqrt(xi) * exp(-abs(y5));
        xiMaxLog.fill( log10(xi1) );
        xiMinLog.fill( log10(xi2) );
        expGap3 = (proc[3].p() + proc[5].p()).m2Calc() / proc[5].m2();
        expGap4 = (proc[4].p() + proc[5].p()).m2Calc() / proc[5].m2();
        yGap.fill( log( max(expGap3, expGap4) ) );
        yGapMin.fill( log( min(expGap3, expGap4) ) );
      }

    // End of event loop.
    }

    // Statistics. SD cross section with xi < 0.05.
    pythia.stat();
    if (hasSD) {
      double sigSD  = pythia.info.sigmaGen(103) + pythia.info.sigmaGen(104);
      double nSD    = pythia.info.nAccepted(103) + pythia.info.nAccepted(104);
      double sigSD005 = sigSD * nxi005 / nSD;
      cout << "\n sigma_SD(both sides) = " << fixed << setprecision(3)
           << sigSD << "\n whereof sigma_SD(xi < 0.05) = "<< sigSD005 << endl;
    }

    // Normalize and print histograms.
    double sigNorm = pythia.info.sigmaGen() / nEvent;
    tLog      *= 10. * sigNorm;
    xiMaxLog  *= 10. * sigNorm;
    xiMinLog  *= 10. * sigNorm;
    xiProdLog *= 10. * sigNorm;
    yGap      *=  4. * sigNorm;
    yGapMin   *=  4. * sigNorm;
    tAbsAcc /= tAbsAll;
    lgxiAcc /= lgxiAll;
    cout << tLog << xiMaxLog << xiMinLog << xiProdLog << yGap << yGapMin;
    cout << nTries << lgxiAll << lgxiAcc << tAbsAll << tAbsAcc;

    // Reset histograms at end of loop over fluxes.
    tLog.null();
    xiMaxLog.null();
    xiMinLog.null();
    xiProdLog.null();
    yGap.null();
    yGapMin.null();
    nTries.null();
    lgxiAll.null();
    lgxiAcc.null();
    tAbsAll.null();
    tAbsAcc.null();
  }

  // Done.
  return 0;
}
