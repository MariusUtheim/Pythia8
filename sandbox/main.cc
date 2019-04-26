#include <chrono>
#include <iostream>

#include "Pythia8/ParticleWidths.h"
#include "Pythia8/LowEnergyHadHad.h"
#include "Pythia8/LowEnergySigma.h"

#include "PDGSigma.h"

void singleEvent()
{
  Pythia pythia;
  pythia.readFile("mymain.cmnd");
  pythia.init();

  pythia.next();
}

void manyEvents()
{
  Pythia pythia;
  pythia.readFile("mymain.cmnd");
  pythia.init();

  int nEvent = pythia.mode("Main:numberOfEvents");
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
    pythia.next();
  
}

void fingerprint()
{
  Pythia pythia("../share/Pythia8/xmldoc", false);
  pythia.readFile("mymain.cmnd");
  pythia.readString("Print:quiet = on");
  pythia.init();

  pythia.next();

  cout << pythia.event[pythia.event.size() - 1].p();
}



constexpr double mN = 0.938;
constexpr double mpi = 0.135;
double pLab(double eCM) {
  if (eCM < 2 * mN) 
    return 0.;
  return sqrt(eCM * eCM * (pow2(eCM / (2. * mN)) - 1));
}

double pLabInv(double mBeam, double mTarget, double p) {
  return sqrt(mBeam * mBeam + mTarget * mTarget 
            + mTarget * sqrt(mBeam * mBeam + p * p));
}


int main(int argc, const char *argv[]) {

  Pythia pythia("../share/Pythia8/xmldoc", false);
  pythia.readFile("mymain.cmnd");
  pythia.init();

  LowEnergyResonance resonance;
  resonance.initPtr(&pythia.rndm, &pythia.particleData);
  if (!resonance.init("../share/Pythia8/xmldoc/ParticleWidths.xml"))
    cout << "Failed to initialize " << endl;

  // Pythia cross sections
  LowEnergySigma sigma;
  sigma.initPtr(&pythia.info, &pythia.particleData, &resonance);

  double eLeft = 0., eRight = 2.;
  double mN = pythia.particleData.m0(2212), mpi = pythia.particleData.m0(111);
  double dm = 0;//2 * mN + mpi;

  Hist sigmaTotal = Hist::plotFunc(
    [&](double e0) { return sigma.sigmaTotal(111, 111, e0); }, 
    "total", 300, eLeft, eRight);

  Hist sigmaEl = Hist::plotFunc(
    [&](double e0) { return sigma.sigmaPartial(111, 111, e0, 2); }, 
    "elastic", 100, eLeft, eRight);

  Hist sigmaDiff = Hist::plotFunc(
    [&](double e0) { return sigma.sigmaPartial(111, 111, e0, 1); }, 
    "diffractive", 100, eLeft, eRight);  

  Hist sigmaRes = Hist::plotFunc(
    [&](double e0) { return sigma.sigmaPartial(111, 111, e0, 7); }, 
    "resonant", 100, eLeft, eRight);

  Hist sigmarho = Hist::plotFunc(
    [&](double e0) { return sigma.sigmaPartial(111, 111, e0, 113); }, 
    "rho", 300, eLeft, eRight);

  Hist sigmaf0 = Hist::plotFunc(
    [&](double e0) { return resonance.getPartialResonanceSigma(111, 111, 9000221, e0); }, 
    "f_0(500)", 300, eLeft, eRight);

  Hist sigmaf2 = Hist::plotFunc(
    [&](double e0) { return sigma.sigmaPartial(111, 111, e0, 225); }, 
    "f_2", 100, eLeft, eRight);

/*
  // PDG cross sections
  Hist sigmaTotalPDG = Hist::plotFunc(
    [&](double e0) { return interpol(ppTotal, pLab(e0 + dm)); },
    "Total (PDG)", 100, eLeft, eRight);

  Hist sigmaElPDG = Hist::plotFunc(
    [&](double e0) { return interpol(ppElastic, pLab(e0 + dm)); },
    "Elastic (PDG)", 100, eLeft, eRight
  );
*/
  
  
  HistPlot plt("myplot");
  plt.frame("myplot");
  plt.add(sigmaTotal, "-");
  plt.add(sigmaEl, "-");
  plt.add(sigmaDiff, "-");
  plt.add(sigmaRes, "-");
  plt.add(sigmarho, "-");
  plt.add(sigmaf0, "-");
  plt.add(sigmaf2, "-");

  //plt.add(sigmaEl, "-");
  //plt.add(sigmaDiff, "-");
  //plt.add(sigmaTotalPDG, "h");
  //plt.add(sigmaElPDG, "h");
  plt.plot();

}
