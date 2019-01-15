#include <chrono>
#include <iostream>
#include "tests.h"
#include "Pythia8/MassDependentWidth.h"
#include "Pythia8/LowEnergyController.h"

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

static Hist plotFunction(string title, double xMin, double xMax, std::function<double(double)> f) {
  double dx = (xMax - xMin) / 100;
  Hist result(title, 100, xMin, xMax);
  for (int i = 0; i < 100; ++i) {
    double x = xMin + (i + 0.5) * dx;
    result.fill(x, f(x));
  }
  cout << endl << endl;
  return result;
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


  LowEnergyData lowEnergyData = pythia.lowEnergyData;

  if (!lowEnergyData.sanityCheck())
    return 1;

  lowEnergyData.print();

//  resonanceData.initPtr(&pythia.particleData, &pythia.massDependentWidths);
//  ifstream stream("ResonanceData.xml");
//  resonanceData.readXML(stream);

//  int c = 0;
//  Hist sigma("x", 101, -0.5, 100.5);
//  for (int x = 0; x <= 100; x++)
//    cout << (sigma.fill(x, x), x) << " " ;

  Hist sigma = plotFunction("Sigma total", 1., 3., 
    [&pythia](double x) { return pythia.lowEnergyController.getTotalSigma(2212, 211, x); }
  );
  cout << sigma;
}