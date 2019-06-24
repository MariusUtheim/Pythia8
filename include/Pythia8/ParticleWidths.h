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

private:

  typedef pair<int, int> keyType;

  class Channel {
  public:
    Channel(Interpolator brIn, int idAIn, int idBIn, int lTypeIn) 
      : br(brIn), idA(idAIn), idB(idBIn), lType(lTypeIn) {}
    Interpolator br;
    int idA, idB;
    int lType;
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
    map<keyType, Channel> decayChannels;
  };

  Info* infoPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  map<int, Entry> entries;

  keyType getKey(int& idR, int idA, int idB) const;

  bool _pickMass1(int idRes, double eCM, double mB, int lType, double& mAOut);

  bool _pickMass2(int id1, int id2, double eCM, int lType,
    double& m1Out, double& m2Out);

  bool _pickMasses(int idA, int idB, double eCM,
    double& mAOut, double& mBOut, int lType);

};

}

#endif