#include <chrono>
#include <iostream>
#include "tests.h"
#include "Pythia8/MassDependentWidth.h"
#include "Pythia8/SigmaResonance.h"
#include "Pythia8/ResonanceData.h"

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

class Test {
public:

  static int count;

  int val;
  void speak() { cout << val << endl; }

  Test() : val(++count) { cout << "()"; speak(); }
  Test(int val) : val(val) { cout << "++"; speak(); }

  Test(const Test& other) : val(other.val) { cout << "!!"; speak(); }
  Test(Test&& other) : val(other.val) { cout << "&&"; speak(); }
};

int Test::count = 10;

int main(int argc, const char *argv[]) {

  Pythia pythia(false);

  ParticleData& particleData = pythia.particleData;

  ResonanceData resonanceData;
  resonanceData.initPtr(&pythia.particleData);
  ifstream stream("ResonanceData.xml");
  resonanceData.readXML(stream);

  if (!resonanceData.sanityCheck())
    return 1;

  resonanceData.print();

  int inA = -321, inB = 2212;
  
  //cout << resonanceData.getIso3(13112) << endl;

  //cout << resonanceData.getResonanceSigma(inA, inB, 1.8) << endl;
  
  
  Hist sigmaRes("sigma(resonance)", 200, 1., 2.5);
  double dx = (2.5 - 1.) / 200;

  for (double eCM = 1. + dx / 2; eCM < 2.5; eCM += dx) 
    sigmaRes.fill(eCM, resonanceData.getResonanceSigma(inA, inB, eCM));

  cout << sigmaRes;

  HistPlot plt("myplot");
  plt.plotFrame("outplot", sigmaRes, "Sigma", "$\\sqrt{s}$", "$\\sigma$", "-");
  

  /*
  Hist sigmaAnn("Sigma(annihilation)", 100, 2., 20.);
  double dx = (20. - 2.) / 100;

  for (double eCM = 2. + dx / 2; eCM < 20.; eCM += dx) 
    sigmaAnn.fill(eCM, resonanceData.getAnnihilationSigma(inA, inB, eCM));

  cout << sigmaAnn;
  */
}