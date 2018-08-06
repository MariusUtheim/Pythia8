// ParticleData.h is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the classes containing particle data.
// DecayChannel contains info on a single decay channel.
// ParticleDataEntry contains info on a single particle species.
// ParticleData collects info on all particles as a map.

#ifndef Pythia8_CrossSectionData_H
#define Pythia8_CrossSectionData_H

#include "Pythia8/Basics.h"
#include "Pythia8/Info.h"
#include "Pythia8/PythiaStdlib.h"
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/Settings.h"
#include "Pythia8/StandardModel.h"

namespace Pythia8 {

//==========================================================================

// Forward reference to some classes.
class InteractionChannel;
class CrossSectionDataEntry;
class CrossSectionData;

//==========================================================================

// This class holds info on a particular interaction channel.

class InteractionChannel {

public:

	InteractionChannel(int onModeIn = 0, double bRatioIn = 0., int meModeIn = 0,
                     vector<int> productsIn = vector<int>())
    : onModeSave(onModeIn), bRatioSave(bRatioIn), meModeSave(meModeIn),
      productsSave(productsIn)
		  { }

	// @TODO: Implement copy constructor

	// Member functions for output.
	int 	 onMode()       const { return onModeSave; } 	
	double bRatio()       const { return bRatioSave; }
  int    meMode()       const { return meModeSave; }
  int    nProducts()    const { return productsSave.size(); }
  int    product(int i) const { return (i >= 0 && i < (int)productsSave.size()) 
                                       ? productsSave[i] : 0; }
  
  const vector<int>& products() { return productsSave; }
  
  bool isActive(int idSgn) const { 
    return (onModeSave == 1 || (onModeSave == 2 && idSgn > 0)
                            || (onModeSave == 3 && idSgn < 0));
  }

private:

  ///Comments
	int onModeSave;
	double bRatioSave;
  int meModeSave;
	vector<int> productsSave; // @TODO Should it be const or otherwise encapsulated?

};

//==========================================================================

// This class holds info on cross-sections for a particular particle pair.

class CrossSectionDataEntry {

public:

	CrossSectionDataEntry(int idAIn = 0, int idBIn = 0, double sigmaIn = 0.)
		: idASave(idAIn), idBSave(idBIn), sigmaSave(sigmaIn) {}

	// @TODO: Implement copy constructor

	void initPtr(CrossSectionData* crossSectionDataPtrIn) {
    crossSectionDataPtr = crossSectionDataPtrIn; }
	
	// Accessors
  int 	 idA() 	 const { return idASave; }
	int 	 idB() 	 const { return idBSave; }
	double sigma() const { return sigmaSave; }

  void idA(int idAIn)        { idASave = idAIn; }
  void idB(int idBIn)        { idBSave = idBIn; }
  void sigma(double sigmaIn) { sigmaSave = sigmaIn; }

	// Reset to empty interaction table.
  void clearChannels() { channels.resize(0); }

  // Add a decay channel to the decay table.
  void addChannel(int onMode = 0, double bRatio = 0., int meMode = 0,
                  vector<int> products = vector<int>()) {
    channels.push_back(InteractionChannel(onMode, bRatio, meMode, products)); }

  // Decay table size.
  int sizeChannels() const { return channels.size(); }

  // Gain access to a channel in the decay table.
  InteractionChannel& channel(int i) { return channels[i]; }
  const InteractionChannel& channel(int i) const { return channels[i]; }

  // Rescale sum of branching ratios to unity.
  void rescaleBR(double newSumBR = 1.);

  // @TODO: This is inconsistent with ParticleData.h convention
  bool hasChanged() const { return hasChangedSave; }
  void hasChanged(bool hasChangedIn) { hasChangedSave = hasChangedIn; }
  
  // Random choice of interaction channel according to branching ratios.
  InteractionChannel& pickChannel();

  // @TODO: Somethign with resonance widths?
  ResonanceWidths* resonancePtr() { return resonancePtrSave; }
  void resonancePtr(ResonanceWidths* resonancePtrIn);
  

private:

	int idASave, idBSave;
	double sigmaSave;
	vector<InteractionChannel> channels;

  // Summed branching ratio of currently open channels.
  double currentBRSum;

  // Pointer to ResonanceWidths object; only used for some particles.
  ResonanceWidths* resonancePtrSave;

  // Pointer to the full cross-section data table.
  CrossSectionData* crossSectionDataPtr;

  bool hasChangedSave;
};

//==========================================================================

// This class holds a map of all CrossSectionDataEntries.

class CrossSectionData {

public:

	CrossSectionData() : infoPtr(0), settingsPtr(0), rndmPtr(0), 
		couplingsPtr(0), isInit(false) {}

	// @TODO: Initialize copy operator

	// Initialize pointers.
  void initPtr(Info* infoPtrIn, Settings* settingsPtrIn, Rndm* rndmPtrIn,
    Couplings* couplingsPtrIn) {infoPtr = infoPtrIn;
    settingsPtr = settingsPtrIn; rndmPtr = rndmPtrIn;
    couplingsPtr = couplingsPtrIn;}

	// Read in database from specific file.
  bool init(string startFile = "../share/Pythia8/xmldoc/ParticleData.xml") {
    initCommon(); return readXML(startFile); }

  // Read and process or list whole (or part of) database from/to an XML file.
  bool readXML(string inFile, bool reset = true);
  bool readXML(istream& is, bool reset = true);
  void listXML(string outFile);// @TODO: Implement 

  bool copyXML(const CrossSectionData &crossSectionDataIn);// @TODO: Implement 

  // @TODO: Implement ...FF functions

  // Keep track whether any readings have failed, invalidating run setup.
  bool readingFailed() { return readingFailedSave; }

  // Print out table of whole database, or of only part of it.
  void listAll() { list(false, true); }
  void listChanged(bool changedRes = false) { list(true, changedRes); }
  void list(bool changedOnly = false, bool changedRes = true);

  // Add new entry.
  CrossSectionDataEntry& addCrossSection(int idAIn, int idBIn, double sigmaIn) {
    pair<int, int> ids = pairify(idAIn, idBIn);
    pdt[ids] = CrossSectionDataEntry(idAIn, idBIn, sigmaIn);
    pdt[ids].initPtr(this);
    return pdt[ids];
  }

  // Query existence of an entry.
  bool isCrossSection(int idAIn, int idBIn) const { // @TODO Different function name?
    pair<int, int> ids = pairify(idAIn, idBIn);
    map<pair<int, int>, CrossSectionDataEntry>::const_iterator found // @TODO: auto?
      = pdt.find(ids);
    if (found == pdt.end()) return false;
    else return true;
  }

  // Query existence of an entry and return a const iterator.
  // @TODO: Return a reference instead of pointer
  CrossSectionDataEntry* findCrossSection(int idAIn, int idBIn) {
    pair<int, int> ids = pairify(idAIn, idBIn);
    map<pair<int, int>, CrossSectionDataEntry>::iterator found = pdt.find(ids);
    if (found == pdt.end()) return NULL;
    return &(found->second);
  }

  // Query existence of an entry and return a const iterator.
  // @TODO: Return a reference instead of pointer
  const CrossSectionDataEntry* findCrossSection(int idAIn, int idBIn) const {
    pair<int, int> ids = pairify(idAIn, idBIn);
    map<pair<int, int>, CrossSectionDataEntry>::const_iterator found 
      = pdt.find(ids);
    if (found == pdt.end()) return NULL;
    return &(found->second);
  }

  // Accessors
  double sigma(int idA, int idB) const {
    const CrossSectionDataEntry* entry = findCrossSection(idA, idB);
    if (entry) return entry->sigma(); 
    else return 0.0;
  }

  void sigma(int idA, int idB, double sigmaIn) {
    CrossSectionDataEntry *entry = findCrossSection(idA, idB);
    if (entry) entry->sigma(sigmaIn);
  }

private:

  // The individual particle need access to the full database.
  friend class CrossSectionDataEntry;

	// Pointer to various information on the generation.
  Info*     infoPtr;

  // Pointer to the settings database.
  Settings* settingsPtr;

  // Pointer to the random number generator.
  Rndm*     rndmPtr;

  // Pointer to Standard Model couplings.
  Couplings*   couplingsPtr;

  // Create a pair of indices, such that the first index is smallest @TODO: Verify
  pair<int, int> pairify(int idAIn, int idBIn) const {
    return idAIn <= idBIn ? pair<int, int>(idAIn, idBIn) 
                          : pair<int, int>(idBIn, idAIn);
  }

  // All particle data stored in a map.
  map<pair<int, int>, CrossSectionDataEntry> pdt; // @TODO: name

  // Pointer to current cross-section (e.g. when reading decay channels).
  CrossSectionDataEntry* crossSectionPtr;

  // Flag that initialization has been performed; whether any failures.
  bool   isInit, readingFailedSave;

  // Method for common setting of particle-specific info.
  void   initCommon();

  // @TODO: Pure functions can be made static 

  // Useful functions for string handling.
  bool   boolString(string tag) { string tagLow = toLower(tag);
    return ( tagLow == "true" || tagLow == "1" || tagLow == "on"
    || tagLow == "yes" || tagLow == "ok" ); }

  // Extract XML value following XML attribute.
  string attributeValue(string line, string attribute);
  bool   boolAttributeValue(string line, string attribute);
  int    intAttributeValue(string line, string attribute);
  double doubleAttributeValue(string line, string attribute);

  // Auxiliary functions to readXML() and copyXML().
  bool loadXML(istream& is, bool reset = true);
};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_CrossSectionData_H

// @TODO: Implement UserHooks
