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

  Pythia pythia("../share/Pythia8/xmldoc", false);

  ParticleData& particleData = pythia.particleData;

  ResonanceData resonanceData;
  resonanceData.initPtr(&pythia.particleData);
  ifstream stream("ResonanceData.xml");
  resonanceData.readXML(stream);

  if (!resonanceData.sanityCheck())
    return 1;

  resonanceData.print();

  int inA = -2212, inB = 2212;
  
  Hist sigmaRes("sigma(resonance)", 200, 0., 2.5);
  double dx = (2.5 - 0.) / 200;

  for (double eCM = 0. + dx / 2; eCM < 2.5; eCM += dx) 
    sigmaRes.fill(eCM, resonanceData.getTotalSigma(inA, inB, eCM));

  cout << sigmaRes;

  HistPlot plt("myplot");
  plt.plotFrame("outplot", sigmaRes, "Sigma", "$\\sqrt{s}$", "$\\sigma$", "-");

}