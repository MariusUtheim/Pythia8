#include <chrono>
#include <iostream>
#include "tests.h"
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

static Hist plotFunction(string title, double xMin, double xMax, std::function<double(double)> f) {
  int nBins = 101;
  double dx = (xMax - xMin) / (nBins - 1);
  Hist result(title, nBins, xMin - 0.5 * dx, xMax + 0.5 * dx);
  for (double x = xMin; x <= xMax + 0.5 * dx; x += dx) {
    cout << x << " --> " << f(x) << endl;
    result.fill(x, f(x));
  }
  cout << endl << endl;
  return result;
}
static Interpolator ppiDiffData(1.9, 3.19642, {
    0., 0.597966, 1.6208, 2.64363, 3.66647, 4.6893, 5.71213, 6.67697,
    7.63029, 8.79865, 10.016, 11.3516, 12.4208, 13.3126, 13.984, 14.5885,
    15.0755, 15.4614, 15.7966, 16.0741, 16.3701, 16.6786, 16.9763,
    17.2395, 17.4756, 17.7065, 17.9797, 18.2733, 18.5514, 18.7887,
    18.9995, 19.1861, 19.3658, 19.4813, 19.5585, 19.6249, 19.6775,
    19.7497, 19.858, 20.0074, 20.1426, 20.2455, 20.3198, 20.3758,
    20.4241, 20.4542, 20.4502, 20.4389, 20.4559, 20.5093, 20.5774,
    20.6341
  });

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

  double eLeft = 0.1, eRight = 5.;
  double mN = pythia.particleData.m0(2212), mL = pythia.particleData.m0(3122);

  Hist sigmaTotal = Hist::plotFunc(
    [&](double p) { return sigma.sigmaTotal(2212, 3122, pLabInv(mL, mN, p)); }, 
    "total", 200, eLeft, eRight);

  Hist sigmaEl = Hist::plotFunc(
    [&](double p) { return sigma.sigmaPartial(2212, 3122, pLabInv(mL, mN, p), 2); }, 
    "elastic", 100, eLeft, eRight);

  Hist sigmaDiff = Hist::plotFunc(
    [&](double p) { return sigma.sigmaPartial(2212, 3122, pLabInv(mL, mN, p), 1); }, 
    "diffractive", 100, eLeft, eRight);


  // PDG cross sections
  Hist sigmaTotalPDG = Hist::plotFunc(
    [&](double p) { return interpol(LambdapTotal, p); },
    "Total (PDG)", 100, eLeft, eRight);

  Hist sigmaElPDG = Hist::plotFunc(
    [&](double p) { return interpol(LambdapElastic, p); },
    "Elastic (PDG)", 100, eLeft, eRight
  );
  
  
  
  HistPlot plt("myplot");
  plt.frame("myplot");
  plt.add(sigmaTotal, "-");
  plt.add(sigmaEl, "-");
  plt.add(sigmaDiff, "-");
  plt.add(sigmaTotalPDG, "h");
  plt.add(sigmaElPDG, "h");
  plt.plot();

}