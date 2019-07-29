#ifndef Hadron_Widths_H
#define Hadron_Widths_H

#include "Pythia8/Info.h"
#include "Pythia8/Interpolator.h"
#include "Pythia8/ParticleData.h"

namespace Pythia8 {

class HadronWidths {

public:

  HadronWidths() = default;
  HadronWidths(const HadronWidths&) = delete;
  HadronWidths(HadronWidths&&) = delete;

  bool init(Info* infoPtrIn, Rndm* rndmPtrIn, ParticleData* particleDataPtrIn,
    string path);

  bool readXML(istream& stream);

  bool check();

  // Returns whether the specified particle is handled by HadronWidths
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

  struct DecayChannel {
    Interpolator partialWidth;
    int idA, idB, lType;
  };

  struct ExcitationChannel {
    Interpolator sigma;
    int maskA, maskB;
  };

  struct Entry {
    double mPeak;
    Interpolator width;
    map<keyType, DecayChannel> decayChannels;
  };

  Info* infoPtr;

  Rndm* rndmPtr;

  ParticleData* particleDataPtr;

  map<int, Entry> entries;
  vector<ExcitationChannel> excitationChannels;

  keyType getKey(int& idR, int idA, int idB) const;

  bool getEntryAndChannel(int idR, int idA, int idB,
    Entry& entryOut, DecayChannel& channelOut) const;

  bool pickMass1(int idRes, double eCM, double mB, int lType, double& mAOut);

  bool pickMass2(int id1, int id2, double eCM, int lType,
    double& m1Out, double& m2Out);

  bool _pickMasses(int idA, int idB, double eCM,
    double& mAOut, double& mBOut, int lType);

};

}

#endif