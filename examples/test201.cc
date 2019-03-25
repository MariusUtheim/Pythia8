// test201.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test program of diffractive scattering in various frameworks.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

int main() {

  // Choose pp/ppbar, CM energy and number of events.
  bool   ispbarp     = false;
  double eCM         = 10000.;
  int    nEvent      = 100000;

  // Choose processes to generate.
  bool   hasEl       = false;
  bool   hasSD       = true;
  bool   hasDD       = true;
  bool   hasCD       = true;

  // Some special cases.
  bool   hasCoulomb  = false;
  bool   compareGrid = false;

  // TotEl & Diff: 0 = own; 1 = SaS/DL, 2 = MBR; 3 = ABMST; 4 = RPP.
  int    iModeTotEl  = 3;
  int    iModeDiff   = 3;

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

  // Damping of small rapidity gaps in all diffractive processes.
  bool   ABMSTdampenGap = true;
  double ABMSTygap   = 2.;
  double ABMSTypow   = 5.;

  // Minimal slopes for diffractive scattering.
  bool   ABMSTuseBMin = true;
  double ABMSTbMinSD = 2.;
  double ABMSTbMinDD = 2.;
  double ABMSTbMinCD = 2.;

  // Some frequent combinations.
  double s           = eCM * eCM;
  double mProton     = 0.938272;
  double mu1         = pow2( mProton / eCM);

  // Book histograms.
  Hist tLog("log10(|t|) spectrum", 80, -6., 2.);
  Hist xiMaxLog("diffractive log10(xi_Max) spectrum", 100, -10., 0.);
  Hist xiMinLog("diffractive log10(xi_Min) spectrum", 100, -10., 0.);
  Hist xiProdLog("diffractive log10(xi_1 * xi_2) spectrum", 100, -10., 0.);
  Hist yGap("rapidity gap size", 100, -5., 20.);
  Hist yGapMin("smallest rapidity gap size", 100, 0., 10.);

  // Loop over different fluxes for user-set diffractive scenario.
  // Also other options accessible.
  int nFlux = (iModeDiff == 0) ? 7 : 1;
  for (int iFlux = 1; iFlux <= nFlux; ++iFlux) {

    // Generator.
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
    pythia.settings.mode("SigmaDiffractive:mode", iModeDiff);
    pythia.settings.mode("Diffraction:PomFlux", iFlux);
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

    // Compare cross section values in grid with Christine's programs.
    if (compareGrid) {

      // Extract Pythia databases and feed to local cross section object.
      SigmaTotal sigLocal;
      sigLocal.init( &pythia.info, pythia.settings, &pythia.particleData,
        &pythia.rndm);
      sigLocal.calc( 2212, (ispbarp ? -2212 : 2212), eCM);
      double xiVal[7]  = { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1 };
      double tVal[5]   = { -0.001, -0.01, -0.1, -1., -10. };

       // Single diffraction: run through grid in xi, t.
      if (hasSD) {
	cout << "\n SD cross sections in some phase-space points: \n"
	     << "    xi1          t        sig_SD"  << endl;
	for (int ixi = 0; ixi < 7; ++ixi)
	for (int it = 0; it < 5; ++it) {
	  double xi      = xiVal[ixi];
	  double t       = tVal[it];

	  // Evaluate cross section and print.
          bool tInRange  = sigLocal.tInRange( t/s, 1., mu1, mu1, xi, mu1);
	  double dsigSD  = (tInRange) ? sigLocal.dsigmaSD( xi, t) : 0.;
	  if (dsigSD > 1e-10) cout << scientific << setprecision(3)
            << setw(11) << xi << setw(12) << t << setw(11) << dsigSD << endl;
        }
      }

      // Double diffraction: run through grid in xi1, xi2, t.
      if (hasDD) {
	cout << "\n      DD cross sections in some phase-space points: \n"
	     << "    xi1        xi2          t        sig_DD     sig_SD1"
	     << "    sig_SD2    sig_el      DD(t)    DD(0)*exp" << endl;
	for (int ixi1 = 0; ixi1 < 7; ++ixi1)
	for (int ixi2 = 0; ixi2 <= ixi1; ++ixi2)
	for (int it = 0; it < 5; ++it) {
	  double xi1     = xiVal[ixi1];
	  double xi2     = xiVal[ixi2];
	  double t       = tVal[it];

	  // Evaluate cross section components and print.
	  // Extra printing if minimal exponential fall-off kicks in.
          bool tInRange  = sigLocal.tInRange( t/s, 1., mu1, mu1, xi1, xi2);
	  double dsigDD  = (tInRange) ? sigLocal.dsigmaDD( xi1, xi2, t) : 0.;
	  if (dsigDD > 1e-10) {
	    double dsigSD1 = sigLocal.dsigmaSD( xi1, t);
	    double dsigSD2 = sigLocal.dsigmaSD( xi2, t);
	    double dsigEl  = sigLocal.dsigmaEl( t, false, false);
	    double dsigDD1 = dsigSD1 * dsigSD2 / dsigEl;
	    double dsigDD2 = (sigLocal.dsigmaSD( xi1, 0.)
	      * sigLocal.dsigmaSD( xi2, 0.)
	      / sigLocal.dsigmaEl( 0., false, false)) * exp(ABMSTbMinDD * t);
	    cout << scientific << setprecision(3) << setw(11) << xi1
		 << setw(11) << xi2 << setw(12) << t << setw(11) << dsigDD
		 << setw(11) << dsigSD1 << setw(11) << dsigSD2
		 << setw(11) << dsigEl;
	    if (dsigDD1 > dsigDD2) cout << setw(11) << dsigDD1
		 << setw(11) << dsigDD2;
	    cout << endl;
	  }
	}
      }

    // End of comparison of values in a grid.
    }

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

      // Study t, xiMax and xiMin distributions of elastic/diffractive events.
      int code = pythia.info.code();
      double tAbsL = log10(abs(pythia.info.tHat()));
      tLog.fill( tAbsL);
      double mMax, xiMx, xiMxL, mMin, xiMnL, yGp, xiProd, xiProdL, y5,
             yGp3, yGp4;
      if (code > 102 && code < 106) {
        mMax  = max( proc[3].m(), proc[4].m());
        xiMx  = pow2( mMax / eCM );
        xiMxL = log10( xiMx );
        xiMaxLog.fill( xiMxL );
        if ((code == 103 || code == 104) && xiMx < 0.05) ++nxi005;
      }
      if (code == 105) {
        mMin  = min( proc[3].m(), proc[4].m());
        xiMnL = log10( pow2( mMin / eCM ) );
        xiMinLog.fill( xiMnL);
        yGp   = log( pow2( eCM * mProton / (mMax * mMin) ) );
        yGap.fill( yGp );
      } else if (code == 106) {
        double uAbsL = log10(abs(pythia.info.uHat()));
        tLog.fill( uAbsL);
        xiProd  = pow2( proc[5].m() / eCM);
        xiProdL = log10( xiProd );
        xiProdLog.fill( xiProdL );
        y5   = proc[5].y();
        xiMxL = log10(sqrt(xiProd) * exp( abs(y5)));
        xiMnL = log10(sqrt(xiProd) * exp(-abs(y5)));
        xiMaxLog.fill( xiMxL );
        xiMinLog.fill( xiMnL);
        yGp3 = log( (proc[3].p() + proc[5].p()).m2Calc() / proc[5].m2() );
        yGap.fill( yGp3 );
        yGp4 = log( (proc[4].p() + proc[5].p()).m2Calc() / proc[5].m2() );
        yGap.fill( yGp4 );
        yGapMin.fill( min(yGp3, yGp4) );
      }

    // End of event loop.
    }

    // Final statistics. Normalize and print histograms.
    pythia.stat();
    if (hasSD) {
      double sigSD  = pythia.info.sigmaGen(103) + pythia.info.sigmaGen(104);
      double nSD    = pythia.info.nAccepted(103) + pythia.info.nAccepted(104);
      double sigSD005 = sigSD * nxi005 / nSD;
      cout << "\n sigma_SD(both sides) = " << fixed << setprecision(3)
           << sigSD << "\n sigma_SD( xi < 0.05) = "<< sigSD005 << endl;
    }
    double sigNorm = pythia.info.sigmaGen() / nEvent;
    tLog     *= 10. * sigNorm;
    xiMaxLog *= 10. * sigNorm;
    xiMinLog *= 10. * sigNorm;
    xiProdLog *= 10. * sigNorm;
    yGap     *=  4. * sigNorm;
    yGapMin  *= 10. * sigNorm;
    cout << tLog << xiMaxLog << xiMinLog << xiProdLog << yGap << yGapMin;

    // Reset histograms at end of loop over fluxes.
    tLog.null();
    xiMaxLog.null();
    xiMinLog.null();
    yGap.null();
  }

  // Done.
  return 0;
}
