#ifndef Particle_Widths_H
#define Particle_Widths_H

#include "Pythia8/Info.h"
#include "Pythia8/Interpolator.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {


//TS?? HadronWidths better name?
//Then also ParticleWidths.xml -> HadronWidthData.xml
class ParticleWidths {

public:

  ParticleWidths() = default;
  ParticleWidths(const ParticleWidths&) = delete;
  ParticleWidths(ParticleWidths&&) = delete;

  bool init(Info* infoPtrIn, Rndm* rndmPtrIn, ParticleData* particleDataPtrIn,
    string path);

  bool readXML(istream& stream);

  // Returns whether the specified particle is handled by ParticleWidths
  bool hasData(int id) const {
    auto iter = entries.find(abs(id));
    return iter != entries.end();
  }

  // Get a list of all implemented resonances
  vector<int> getResonances() const;

  // Get the total width of the specified particle at specified mass
  double width(int id, double m) const;

  // Get the partial width of the specified particle and products
  double partialWidth(int id, int prodA, int prodB, double m) const;

  // Get the branching ratio of the specified particle at specified mass
  double branchingRatio(int id, int prodA, int prodB, double m) const;

  // Calculate resonance formation cross section
  double resonanceSigma(int idA, int idB, int idRes, double eCM) const;

  // Pick masses for two particles, taking distributions and phase space
  // into account. Returns whether successful.
  bool pickMasses(int idA, int idB, double eCM, double& mAOut, double& mBOut);

  // Pick a decay channel for the specified particle, together with phase
  // space configuration. Returns whether successful.
  bool pickDecay(int idDec, double m, int& id1Out, int& id2Out,
    double& m1Out, double& m2Out);

  bool pickExcitation(int idA, int idB, double eCM, 
    int& idCOut, double& mCOut, int& idDOut, double& mDOut);

private:

  typedef pair<int, int> keyType;

  class DecayChannel {
  public:
    DecayChannel(Interpolator brIn, int idAIn, int idBIn, int lTypeIn) 
      : br(brIn), idA(idAIn), idB(idBIn), lType(lTypeIn) {}
    DecayChannel(const DecayChannel&) = delete; // prohibit copying
    DecayChannel(DecayChannel&&) = default;
    Interpolator br;
    int idA, idB, lType;
  };

  class ExcitationChannel {
  public:
    ExcitationChannel(Interpolator sigmaIn, int maskAIn, int maskBIn) 
      : sigma(sigmaIn), maskA(maskAIn), maskB(maskBIn) {}
    Interpolator sigma;
    int maskA, maskB;
  };

  class Entry {
  public:
    Entry(double m0In, Interpolator widthsIn) : m0(m0In), widths(widthsIn) {}
    Entry(const Entry&) = delete;
    Entry(Entry&&) = default;

    void addProducts(keyType, Interpolator brs, int idA, int idB, int lType);
    double getWidth(keyType key, double eCM) const;
    double getBR(keyType key, double eCM) const; 
    int getlType(keyType key) const;

    double m0;
    Interpolator widths;
    map<keyType, DecayChannel> decayChannels;
  };

  Info* infoPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  map<int, Entry> entries;
  vector<ExcitationChannel> excitationChannels;

  keyType getKey(int& idR, int idA, int idB) const;

  bool _pickMass1(int idRes, double eCM, double mB, int lType, double& mAOut);

  bool _pickMass2(int id1, int id2, double eCM, int lType,
    double& m1Out, double& m2Out);

  bool _pickMasses(int idA, int idB, double eCM,
    double& mAOut, double& mBOut, int lType);

};

}

#endif