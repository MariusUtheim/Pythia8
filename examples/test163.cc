// Written by Neelima Sehgal 2012
// program outputs the gamma ray spectrum looping though dark matter masses

//#include "Pythia.h"
#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// A derived class for (e+ e- ->) GenericResonance -> various final states.

class Sigma1GenRes : public Sigma1Process {

public:

  // Constructor.
  Sigma1GenRes() {}

  // Evaluate sigmaHat(sHat): dummy unit cross section.
  virtual double sigmaHat() {return 1.;}

  // Select flavour. No colour or anticolour.
  virtual void setIdColAcol() {setId( -11, 11, 999999);
    setColAcol( 0, 0, 0, 0, 0, 0);}

  // Info on the subprocess.
  virtual string name()    const {return "GenericResonance";}
  virtual int    code()    const {return 9001;}
  virtual string inFlux()  const {return "ffbarSame";}

};

//==========================================================================

int main() {

  // Generator.
  Pythia pythia;

  // A class to generate the fictitious resonance initial state.
  SigmaProcess* sigma1GenRes = new Sigma1GenRes();

  // Hand pointer to Pythia.
  pythia.setSigmaPtr( sigma1GenRes);

  // Read in commands from external file.
  pythia.readFile("test163.cmnd");

  // Extract settings to be used in the main program.
  int  nEvent  = pythia.mode("Main:numberOfEvents");
  int  nList   = 1;

  // Histogram photon spectra.
  Hist eGamma("log10(energy spectrum of photons)", 100, -8., 2.);
  Hist eGammaX("x-weighted kog10(energy) photons", 100, -8., 2.);

  // Set CM energy and initialize.
  double DMenergy = 50.;
  double CME = 2. * DMenergy;
  pythia.settings.forceParm("Beams:eCM", CME);
  pythia.particleData.m0(999999, CME);
  pythia.init();
  int nLowGam = 0;

  // Loop over event generation.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < nList) {pythia.info.list(); pythia.event.list();}

    // Fill histogram.
    for (int i = 0; i < pythia.event.size(); ++i){
      if (pythia.event[i].isFinal() && pythia.event[i].idAbs() == 22) {
        double eI    = pythia.event[i].e();
        double eLog  = log10(eI);
        double xFrac = eI/DMenergy;
        eGamma.fill(eLog);
        eGammaX.fill(eLog, 1./ xFrac);
        if (xFrac < 1e-4) ++nLowGam;
      }
    }

  // End of event loop.
  }

  // Statistics, histogram and done.
  pythia.stat();
  cout << eGamma << eGammaX;
  cout << "\n Number of low-energy photons = " << nLowGam;

  return 0;
}
