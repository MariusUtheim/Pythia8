#ifndef Low_Energy_Controller_H
#define Low_Energy_Controller_H

#include "Basics.h"
#include "Event.h"
#include "Info.h"
#include "ParticleData.h"

namespace Pythia8 {

//==========================================================================

// LowEnergyHadHad: described the low-energy collision between two hadrons.

class LowEnergyHadHad {

public:

  // Constructor. Still to be expanded with further default values.
  LowEnergyHadHad() : infoPtr(), particleDataPtr(), rndmPtr() {}

  // Initialize the class.
  bool init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn);

  // Produce outgoing primary hadrons from collision of incoming pair.
  bool collide( int i1, int i2, int type, Event& event, Vec4 vtx = Vec4() );

  // Event record to handle hadronization.
  Event         leEvent; 

private:

  // Constants: could only be changed in the code itself.
  static const int MAXLOOP;
  static const double MASSREDUCERATE, MDIFFMIN, ALPHAPRIME;

  // Parameters of the generation process.
  double fracEtass, fracEtaPss, xPowMes, xPowBar, xDiqEnhance, sigmaQ;

  // Properties of the current collision. 1 or 2 is two incoming hadrons.
  // "c" or "ac" is colour or anticolour component of hadron.
  bool   isBaryon1, isBaryon2; 
  int    sizeOld, id1, id2, idc1, idac1, idc2, idac2;
  double m1, m2, eCM, sCM, z1, z2, mT1, mT2, mA, mB, 
         mc1, mac1, px1, py1, pTs1, mTsc1, mTsac1, mTc1, mTac1,
         mc2, mac2, px2, py2, pTs2, mTsc2, mTsac2, mTc2, mTac2;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Handle inelastic nondiffractive collision.
  bool nondiff();

  // Handle elastic and diffractive collisions.
  bool eldiff( int type);

  // Handle annihilation collisions.
  bool annihilation();

  // Split up hadron A or B into a colour pair, with masses and pT values.
  bool splitA( double redMpT);
  bool splitB( double redMpT); 

  // Split a hadron inte a colour and an anticolour part.
  pair< int, int> splitFlav( int id);  

  // Choose relative momentum of colour and anticolour constituents in hadron.
  double splitZ( int iq1, int iq2, double mRat1, double mRat2);

  // Estimate lowest possible mass state for flavour combination.
  double mThreshold( int iq1, int iq2);

  // Pick slope b of exp(b * t) for elastic and diffractive events.
  double bSlope( int type);

};


}

#endif