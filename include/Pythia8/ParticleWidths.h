#ifndef Particle_Widths_H
#define Particle_Widths_H

#include "Pythia8/Info.h"
#include "Pythia8/Interpolator.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {

class MDWDecayChannel {
public:

  MDWDecayChannel(Interpolator brIn, int lTypeIn)
    : br(brIn), lType(lTypeIn) {}

  Interpolator br;
  int lType;

};

class ParticleWidthEntry {
public:

  ParticleWidthEntry(double m0In, Interpolator widthsIn)
    : m0(m0In), widths(widthsIn) {}
  ParticleWidthEntry(const ParticleWidthEntry&) = delete;
  ParticleWidthEntry(ParticleWidthEntry&&) = default;

  void addProducts(pair<int, int> prods, Interpolator brs, int lType) {
    decayChannels.emplace(prods, MDWDecayChannel(brs, lType));
  }

  double getWidth(pair<int, int> prods, double eCM) const {
    auto iter = decayChannels.find(prods);
    return (iter != decayChannels.end()) ? iter->second.br(eCM) * widths(eCM) : 0.;
  }

  double getBR(pair<int, int> prods, double eCM) const {
    auto iter = decayChannels.find(prods);
    return (iter != decayChannels.end()) ? iter->second.br(eCM) : 0.;
  }

  int getlType(pair<int, int> prods) const {
    auto iter = decayChannels.find({ abs(prods.first), abs(prods.second) });
    return (iter != decayChannels.end()) ? iter->second.lType : 0;
  }

  double m0;
  Interpolator widths;
  map<pair<int, int>, MDWDecayChannel> decayChannels;
};

//TS?? HadronWidths better name?
//Then also ParticleWidths.xml -> HadronWidthData.xml
class ParticleWidths {

public:

  bool init(Info* infoPtrIn, Rndm* rndmPtrIn, ParticleData* particleDataPtrIn,
    string path) {
    infoPtr = infoPtrIn;
    rndmPtr = rndmPtrIn;
    particleDataPtr = particleDataPtrIn;

    ifstream stream(path);
    if (!stream.is_open()) {
      infoPtr->errorMsg( "Warning in ParticleWidths::init: "
          "unable to open file");
      return false;
    }
    return readXML(stream);
  }

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

  Info* infoPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  map<int, ParticleWidthEntry> entries;

  pair<int, int> getKey(int& idR, int idA, int idB) const;

  bool _pickMass1(int idRes, double eCM, double mB, int lType, double& mAOut);

  bool _pickMass2(int id1, int id2, double eCM, int lType,
    double& m1Out, double& m2Out);

  bool _pickMasses(int idA, int idB, double eCM,
    double& mAOut, double& mBOut, int lType);

};

}

#endif