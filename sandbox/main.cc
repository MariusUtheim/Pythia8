#include <chrono>
#include <iostream>

#include "Pythia8/ParticleWidths.h"
#include "Pythia8/LowEnergyHadHad.h"
#include "Pythia8/LowEnergySigma.h"

#include "PDGSigma.h"



int main(int argc, const char *argv[]) {

  Pythia pythia("../share/Pythia8/xmldoc", false);
  pythia.readFile("mymain.cmnd");

  pythia.init();

  LowEnergyHadHad leHadHad;
  leHadHad.init(&pythia.info, pythia.settings, &pythia.particleData, 
    &pythia.rndm);
  

  double eCM = 1.9;
  double sCM = eCM * eCM;

  int idA = 2212, idB = 211;
  double mA = pythia.particleData.m0(idA);
  double mB = pythia.particleData.m0(idB);

  double pcms = sqrt((sCM - pow2(mA + mB)) * (sCM - pow2(mA - mB))) / (2 * eCM);

  for (int i = 0; i < 10; i++) {
    pythia.event.clear();
    pythia.event.append(idA, 1, 0,0,0,0, 0,0,
                        0., 0.,  pcms, sqrt(pcms * pcms + mA * mA), mA);
    pythia.event.append(idB, 1, 0,0,0,0, 0,0,
                        0., 0., -pcms, sqrt(pcms * pcms + mB * mB), mB);

    leHadHad.collide(0, 1, 0, pythia.event);

    pythia.event.list();
  }
}
