// test179.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Compare pT spectra in thermal vs ordinary hadronization models.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  // Number of events.
  int nEvent = 1000;

  // Histogram multiplicities and pT spectra.
  Hist nChN("charged multiplicity in normal model", 100, 1., 401.);
  Hist nChT("charged multiplicity in thermal model", 100, 1., 401.);
  Hist nChM("charged multiplicity in mT2 model", 100, 1., 401.);
  Hist pTpiN("pT for pi+- in normal model", 100, 0., 5.);
  Hist pTpiT("pT for pi+- in thermal model", 100, 0., 5.);
  Hist pTpiM("pT for pi+- in mT2 model", 100, 0., 5.);
  Hist pTKN("pT for K+- in normal model", 100, 0., 5.);
  Hist pTKT("pT for K+- in thermal model", 100, 0., 5.);
  Hist pTKM("pT for K+- in mT2 model", 100, 0., 5.);
  Hist pTpN("pT for p/pbar in normal model", 100, 0., 5.);
  Hist pTpT("pT for p/pbar in thermal model", 100, 0., 5.);
  Hist pTpM("pT for p/pbar in mT2 model", 100, 0., 5.);

  // Loop over nromal and thermal generation.
  for (int iModel = 0; iModel < 3; ++iModel) {

    // Generator. Common process selection.
    Pythia pythia;
    pythia.readString("Beams:eCM = 8000.");
    pythia.readString("SoftQCD:nonDiffractive = on");

    // Model-specific setup. Initialization.
    if (iModel == 1) pythia.readString("StringPT:thermalModel = on");
    if (iModel == 2) pythia.readString("StringPT:mT2suppression = on");
    pythia.init();
    Event& event = pythia.event;

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Find charged multiplcity and fill histograms.
      int nCharged = 0;
      for (int i = 0; i < event.size(); ++i)
        if (event[i].isFinal() && event[i].isCharged()) ++nCharged;
      if      (iModel == 0) nChN.fill( nCharged );
      else if (iModel == 1) nChT.fill( nCharged );
      else                  nChM.fill( nCharged );

      // Find pT spectra of some hadron species.
      for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {
        int idAbs = event[i].idAbs();
        double pT = event[i].pT();
        if (iModel == 0) {
          if (idAbs ==  211) pTpiN.fill( pT);
          if (idAbs ==  321) pTKN.fill( pT);
          if (idAbs == 2212) pTpN.fill( pT);
        } else if (iModel == 1) {
          if (idAbs ==  211) pTpiT.fill( pT);
          if (idAbs ==  321) pTKT.fill( pT);
          if (idAbs == 2212) pTpT.fill( pT);
        } else {
          if (idAbs ==  211) pTpiM.fill( pT);
          if (idAbs ==  321) pTKM.fill( pT);
          if (idAbs == 2212) pTpM.fill( pT);
        }
      }

    // End of event loop. Statistics.
    }
    pythia.stat();

  // End of model loop. Histograms.
  }
  cout << nChN << nChT << nChM << pTpiN << pTpiT << pTpiM
       << pTKN << pTKT << pTKM << pTpN << pTpT << pTpM;

  // Done.
  return 0;
}
