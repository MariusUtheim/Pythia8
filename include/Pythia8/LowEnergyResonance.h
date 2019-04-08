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

  // Form a resonance between two particles, then decay it
  bool collide(int i1, int i2, Event& event, Vec4 origin = Vec4());

  // Get sigma for resonance scattering through the specified resonance
  double getPartialResonanceSigma(int idA, int idB, int idR, double eCM) const;

  // Get sigma for resonance scattering through all possible resonances
  double getResonanceSigma(int idA, int idB, double eCM) const;

//  double getPartialElasticResonanceSigma(int idA, int idB, int idR, double eCM) const;
//
//  double getElasticResonanceSigma(int idA, int idB, double eCM) const;

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

  // Get a list over possible resonances that can be formed by the particles
  vector<int> getResonanceCandidates(int idA, int idB) const;
};

}

#endif