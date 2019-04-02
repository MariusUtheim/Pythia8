#include <chrono>
#include <iostream>
#include "tests.h"
#include "Pythia8/ParticleWidths.h"
#include "Pythia8/LowEnergyHadHad.h"

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
  pythia.readFile("mymain.cmnd");
  pythia.init();

  LowEnergyHadHad lowEnergyHadHad;
  lowEnergyHadHad.init(&pythia.info, pythia.settings, &pythia.particleData, &pythia.rndm);

  cout << "Setting up beams..." << endl;
  pythia.event.reset();
  pythia.event.append(2212, 12, 0, 0, 0, 0, 0, 0, 
                             0, 0, 4, sqrt(16 + pythia.particleData.m0(2212)));
  pythia.event.append(211, 12, 0, 0, 0, 0, 0, 0,
                             0, 0, -4, sqrt(16 + pythia.particleData.m0(211)));

  cout << "Performing collision..." << endl;
  lowEnergyHadHad.collide(1, 2, 7, pythia.event);

  pythia.event.list(false, false);
}