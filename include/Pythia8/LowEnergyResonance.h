#ifndef Low_Energy_Resonance_H
#define Low_Energy_Resonance_H

#include "Pythia8/Event.h"
#include "Pythia8/ParticleWidths.h"

namespace Pythia8 {

//==========================================================================

// This deals with cross sections and scattering through resonances

class LowEnergyResonance {
public:

  //TS?? Could be combined with init.
  void initPtr(Rndm* rndmPtrIn, ParticleData* particleDataPtrIn)
  { rndmPtr = rndmPtrIn; particleDataPtr = particleDataPtrIn; }

  //TS?? Compare with code for parton distributions, e.g.
  //LHAGrid1(int idBeamIn = 2212, string pdfWord = "void",
  //  string xmlPath = "../share/Pythia8/xmldoc/",
  //where xmlPath is set up in Pythia.cc, the Pythia constructor and saved:
  //settings.addWord( "xmlPath", xmlPath); 
  bool init(string path);

  // Get sigma for resonance scattering through the specified resonance
  double getPartialResonanceSigma(int idA, int idB, int idR, double eCM) const;

  // Get sigma for resonance scattering through all possible resonances
  double getResonanceSigma(int idA, int idB, double eCM) const;

  // Get a list over possible resonances that can be formed by the particles
  // This is determined only by conservation of quantum numbers; just beacuse
  // a resonance conserves quantum numbers, it does not mean it can actually
  // decay into the specified products.
  vector<int> getPossibleResonances(int idA, int idB) const;

  int pickResonance(int idA, int idB, double eCM);
  vector<int> pickDecayProducts(int idRes, double eCM);

private:

  Rndm* rndmPtr;

  Info* infoPtr;

  ParticleData* particleDataPtr;

  ParticleWidths particleWidths;

  // @TODO: Make a more intutive system
  // The signature of a particle is the three digit number BQS, where B is 
  // baryon number, Q is charge signature and S is strangeness signature.
  // A resonance can be formed only if it conserves the total signature. 
  // The charge signature of a particle with charge q is given by chargeType if
  // charge is positive and 10 + chargeType if it is negative. This way, charge
  // signature is always positive. Strangeness signature is defined similarly.
  map<int, vector<int>> signatureToParticles;

};

}

#endif