#ifndef Low_Energy_Strings_H
#define Low_Energy_Strings_H

#include "Pythia8/Basics.h"
#include "Pythia8/Event.h"
#include "Pythia8/Info.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/LowEnergyProcess.h"
#include "Pythia8/Settings.h"

namespace Pythia8 {

//==========================================================================

class LowEnergyStrings : public LowEnergyProcess {

public:

  // Constructor.
  LowEnergyStrings() : infoPtr(), particleDataPtr(), rndmPtr() {}

  // Initialize the class.
  bool init(Info* infoPtrIn, Settings* settingsPtrIn,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn);

  // Produce outgoing primary hadrons frrom collision of incoming pair.
  bool collide( int i1, int i2, Event& event);

private:

  // Constants: could only be changed in the code itself.
  static const int MAXLOOP;
  static const double MASSREDUCERATE;

  // Properties of the current collision.
  int    sizeOld;
  double fracetass, fracetaPss, xPowMes, xPowBar, xDiqEnhance, sigmaQ,
         m1, m2, eCM;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Handle inelastic nondiffractive collision.
  bool nondiff( int i1, int i2, Event& event);

  // Split a hadron inte a colour and an anticolour part.
  pair< int, int> splitHad( int id);  

  // Choose relative momentum of colour and anticolour constituents in hadron.
  double splitZ( int id1, int id2, double mRat1, double mRat2);

};

}

#endif