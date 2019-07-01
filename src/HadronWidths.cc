
#include "Pythia8/HadronWidths.h"

namespace Pythia8 {

typedef pair<int, int> keyType;

// @TODO Clean up these static functions
// @TODO Go through error messages

static string attributeValue(string line, string attribute) {
  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

static int intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;
}

static double doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;
}

static void completeTag(istream& stream, string& line) {
  while (line.find(">") == string::npos) {
    string addLine;
    if (!getline(stream, addLine)) break;
    line += " " + addLine;
  } 
}

// Gets key for the decay and flips idR if necessary
keyType HadronWidths::getKey(int& idR, int idA, int idB) const {

  if (idR < 0) {
    idR = -idR;
    idA = particleDataPtr->antiId(idA);
    idB = particleDataPtr->antiId(idB);
  }

  if (abs(idA) < abs(idB))
    return { idB, idA };
  else
    return { idA, idB };
}


bool HadronWidths::init(Info* infoPtrIn, Rndm* rndmPtrIn,
  ParticleData* particleDataPtrIn, string path) {
  infoPtr = infoPtrIn;
  rndmPtr = rndmPtrIn;
  particleDataPtr = particleDataPtrIn;

  ifstream stream(path);
  if (!stream.is_open()) {
    infoPtr->errorMsg( "Error in HadronWidths::init: "
        "unable to open file");
    return false;
  }
  return readXML(stream);
}

bool HadronWidths::readXML(istream& stream) {

  string line;

  while (getline(stream, line)) {

    string word1;
    istringstream(line) >> word1;

    if (word1 == "<excitationChannel") {
      completeTag(stream, line);

      int maskA = intAttributeValue(line, "maskA");
      int maskB = intAttributeValue(line, "maskB");

      double left  = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");
      
      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      Interpolator sigmas(left, right, data);
      excitationChannels.push_back(ExcitationChannel{ sigmas, maskA, maskB });
    }
    else if (word1 == "<width") {
      completeTag(stream, line);

      int id    = intAttributeValue(line, "id");

      if (entries.find(id) != entries.end()) {
        infoPtr->errorMsg( "Error in HadronWidths::readXML: "
          "resonance is defined more than once",
          std::to_string(id));
        continue;
      }

      double left  = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");
      double m0    = doubleAttributeValue(line, "m0");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      entries.emplace(id, Entry{m0, Interpolator(left, right, data), {}});
    }
    // @TODO rename br to partial width
    else if (word1 == "<br") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");
      auto entryIter = entries.find(id);
      if (entryIter == entries.end()) {
        infoPtr->errorMsg( "Error in HadronWidths::readXML: "
          "got br for a particle with undefined or ill-defined width",
          std::to_string(id));
        continue;
      }

      int lType = intAttributeValue(line, "lType");
      
      istringstream productStr(attributeValue(line, "products"));
      int prod1, prod2;
      productStr >> prod1;
      productStr >> prod2;

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      keyType key = getKey(id, prod1, prod2);
      
      auto& entry = entryIter->second;
      Interpolator br(entry.widths.left(), entry.widths.right(), data);
      entry.decayChannels.emplace(key, DecayChannel{br, prod1, prod2, lType});
    }
  }

  return true;
}

bool HadronWidths::check() {

  // Check that all excitations make sense
  for (auto excitationChannel : excitationChannels) {
    // Check that ids actually correspond to particles
    for (int mask : { excitationChannel.maskA, excitationChannel.maskB })
    for (int id : { mask + 2210, mask + 2110 })
    if (!particleDataPtr->isParticle(id)) {
      infoPtr->errorMsg("Error in HadronWidths::check: "
        "excitation is not a particle", std::to_string(id));
      return false;
    }
  }

  // Check that all resonance entries make sense
  for (auto entryPair : entries) {
    int id = entryPair.first;
    auto& entry = entryPair.second;
    
    // Check that entry id actually corresponds to a particle
    if (!particleDataPtr->isParticle(id)) {
      infoPtr->errorMsg("Error in HadronWidths::check: "
        "resonance is not a particle", std::to_string(id));
      return false;
    }

    // Check that entry id is positive (antiparticles are handled by symmetry)
    if (id < 0) {
      infoPtr->errorMsg("Error in HadronWidths::check: "
        "resonance is an antiparticle", std::to_string(id));
      return false;
    }

    // Check that entry id is positive (antiparticles are handled by symmetry)
    if (!particleDataPtr->isHadron(id)) {
      infoPtr->errorMsg("Error in HadronWidths::check: "
        "resonance is not a hadron", std::to_string(id));
      return false;
    }

    // Check that entry has no charm or bottom quarks
    if (abs(particleDataPtr->heaviestQuark(id)) > 3) {
      infoPtr->errorMsg("Error in HadronWidths::check: "
        "resonance contains a charm or bottom quark", std::to_string(id));
      return false;
    }

    // Check that all decay channels make sense
    for (auto channelPair : entry.decayChannels) {
      auto& channel = channelPair.second;
      int idA = channel.idA, idB = channel.idB;
      string channelStr = std::to_string(id) + " --> " 
          + std::to_string(idA) + " + " + std::to_string(idB);

      // Check that decay product ids actually correspond to particles
      for (int idProd : { idA, idB }) 
      if (!particleDataPtr->isParticle(idProd)) {
        infoPtr->errorMsg("Error in HadronWidths::check: "
          "decay product is not a particle", std::to_string(idProd));
        return false;
      }
      // Check that lType makes sense
      if (channel.lType == 0) {
        infoPtr->errorMsg("Error in HadronWidths::check: "
          "decay channel does not specify an lType", channelStr);
        return false;
      }

      // Check that decay conserves charge (for debugging purposes)
      if (particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB)
        != particleDataPtr->chargeType(id)) {
        infoPtr->errorMsg("Error in HadronWidths::check: "
          "decay does not conserve charge", channelStr);
        return false;
      }
    }
  }

  // @TODO Check that all particles that useMassDependentWidth have data


  return true;

  // @TODO call check() whenever HadronWidths has been initialized


  // Checks on decay channels:
  //  -- Products are particles exists
  //  -- Charge is conserved
  //  -- lType makes sense

  return true;
}

bool HadronWidths::_getEntryAndChannel(int idR, int idA, int idB, 
  const Entry*& entryOut, const DecayChannel*& channelOut) const {

  auto key = getKey(idR, idA, idB);

  auto entryIter = entries.find(idR);
  if (entryIter == entries.end())
    return false;

  auto channelIter = entryIter->second.decayChannels.find(key);
  if (channelIter == entryIter->second.decayChannels.end())
    return false;
  
  entryOut   = &entryIter->second;
  channelOut = &channelIter->second;
  return true;
}


vector<int> HadronWidths::getResonances() const {
  vector<int> resonances;
  for (auto& p : entries)
    resonances.push_back(p.first);
  return resonances;
}

double HadronWidths::width(int id, double eCM) const {
  auto iter = entries.find(abs(id));
  return (iter != entries.end()) ? iter->second.widths(eCM) : 0.;
}

double HadronWidths::branchingRatio(int idR, int idA, int idB, double m) const {
  const Entry* entry;
  const DecayChannel* channel;
  return _getEntryAndChannel(idR, idA, idB, entry, channel)
       ? channel->br(m) : 0.;
}

double HadronWidths::partialWidth(int idR, int idA, int idB, double m) const {
  const Entry* entry;
  const DecayChannel* channel;
  return _getEntryAndChannel(idR, idA, idB, entry, channel)
       ? entry->widths(m) * channel->br(m) : 0.;
}


bool HadronWidths::pickMasses(int idA, int idB, double eCM, double& mAOut, double& mBOut) {
  return _pickMasses(idA, idB, eCM, mAOut, mBOut, 1);
}

bool HadronWidths::_pickMasses(int idA, int idB, double eCM,
  double& mAOut, double& mBOut, int lType) {

  bool isResA = hasData(idA), isResB = hasData(idB);

  // If neither particle has mass-dependent width
  if (isResA && isResB) {
    return _pickMass2(idA, idB, eCM, lType, mAOut, mBOut);
  }
  else if (isResA) {
    double mA, mB = particleDataPtr->m0(idB);
    if (!_pickMass1(idA, eCM, mB, lType, mA))
      return false;
    mAOut = mA; mBOut = mB;
    return true;
  }
  else if (isResB) {
    double mB, mA = particleDataPtr->m0(idA);
    if (!_pickMass1(idB, eCM, mA, lType, mB))
      return false;
    mAOut = mA; mBOut = mB;
    return true;
  }
  else {
    mAOut = particleDataPtr->m0(idA);
    mBOut = particleDataPtr->m0(idB);
    return true;
  }
}

static double pCMS(double eCM, double mA, double mB) {
  double sCM = eCM * eCM;
  return sqrt((sCM - pow2(mA + mB)) * (sCM - pow2(mA - mB))) / (2. * eCM);
}

static double breitWigner(double gamma, double dm) {
  return 1. / (2. * M_PI) * gamma / (dm * dm + 0.25 * gamma * gamma);
}

static constexpr double MAX_LOOPS = 1000;

bool HadronWidths::_pickMass1(int idRes, double eCM, double mB, int lType,
  double& mAOut) {

  // Ensure resonance is positive - the mass distribution doesn't change
  idRes = abs(idRes);

  // Get width entry
  auto iter = entries.find(idRes);
  if (iter == entries.end()) {
    infoPtr->errorMsg("Error in HadronWidths::pickMass: "
      "mass distribution for particle is not defined",
      std::to_string(idRes));
    return false;
  }
  Entry& entry = iter->second;

  // @TODO: Maybe an mPeak that is different from m0 will be more efficient
  double mMin = entry.widths.left(), 
         mMax = min(entry.widths.right(), eCM - mB),
         mPeak = particleDataPtr->m0(idRes),
         m0 = particleDataPtr->m0(idRes);

  // This can happen due to interpolation imprecision if eCM - mB is near mMin
  if (mMax < mMin)
    return false;

  if (mMin > m0 || mMin > mPeak) {
    infoPtr->errorMsg("Error in HadronWidths::pickMass: "
      "mass ordering did not make sense. This indicates an error with the configuration files");
    return false;
  }

  double gamma = entry.widths(mPeak);

  // Scale factor is chosen based on what empirically gives sufficient
  // accuracy and a high efficiency.
  double scale = 1.2 * pow(0.5 * pCMS(eCM, mMin, mB) + 0.5 * pCMS(eCM, mPeak, mB), lType);

  // Get total probability below/above peak
  double cdfLow, cdfHigh;
  if (mPeak < mMax) {
    cdfLow = -1./M_PI * atan((mMin - mPeak) / (0.5 * gamma));
    cdfHigh = 1./M_PI * atan((mMax - mPeak) / gamma);
  }
  else {
    // If max is below peak, i.e. if eCM - mB < m0Res.
    cdfLow = 1./M_PI * atan((mMax - mMin) / (0.5 * gamma));
    cdfHigh = 0.;
  }

  // Relative probabilities of being below/above peak
  vector<double> ps = { cdfLow, 2. * cdfHigh };

  // @TODO: Compare w/ hit-and-miss implementations in other parts of code
  for (int i = 0; i < MAX_LOOPS; ++i) {

    double mCand, envelope;

    if (rndmPtr->pick(ps) == 0) {
      double r = (0.5 - cdfLow) + rndmPtr->flat() * cdfLow;
      mCand = mPeak + 0.5 * gamma * tan(M_PI * (r - 0.5));
      envelope = scale * breitWigner(gamma, mCand - mPeak);
    }
    else {
      double r = 0.5 + rndmPtr->flat() * cdfHigh;
      mCand = mPeak + gamma * tan(M_PI * (r - 0.5));
      envelope = scale * 2. * breitWigner(2. * gamma, mCand - mPeak);
    }
    
    // Rejection step
    double yDistr = pow(pCMS(eCM, mCand, mB), lType) 
                  * breitWigner(entry.widths(mCand), mCand - m0);
  
    if (rndmPtr->flat() * envelope < yDistr) {
      mAOut = mCand;
      return true;
    }
  }

  infoPtr->errorMsg("Warning in HadronWidths::pickMass: "
    "Could not pick mass within prescribed number of iterations. ",
    std::to_string(idRes) + " in (" + std::to_string(mMin) + ", " + std::to_string(mMax) + ")");
  
  mAOut = mMin + rndmPtr->flat() * (mMax - mMin);
  return true;
}

// @TODO: Implement this properly
//        For the first iteration, we just pick id2 on shell
bool HadronWidths::_pickMass2(int id1, int id2, double eCM, int lType,
  double& m1Out, double& m2Out) {
  
  double m2 = particleDataPtr->m0(id2);

  auto iter = entries.find(id2);
  if (iter == entries.end()) return false;

  double m1;
  if (!_pickMass1(id1, eCM, m2, lType, m1))
    return false;

  m1Out = m1; m2Out = m2;
  return true;
}



bool HadronWidths::pickDecay(int idDec, double m, int& idAOut, int& idBOut,
    double& mAOut, double& mBOut) {

  bool isAnti = (idDec < 0);
  if (isAnti)
    idDec = -idDec;

  auto entriesIter = entries.find(idDec);
  if (entriesIter == entries.end()) {
    infoPtr->errorMsg("Error in HadronWidths::pickDecay: "
      "particle not found:", std::to_string(idDec));
    return false;
  }
  const auto& entry = entriesIter->second;

  // Pick decay channel
  vector<pair<int, int>> prodss;
  vector<double> brs;
  bool gotAny = false;
  for (const auto& channel : entry.decayChannels) {
    prodss.push_back(channel.first);
    double br = channel.second.br(m);
    if (br > 0) gotAny = true;
    brs.push_back(br);
  }

  if (!gotAny) {
    infoPtr->errorMsg("Warning in HadronWidths::pickDecay: "
      "no channels have positive branching ratios");
    return false;
  }
  
  auto prods = prodss[rndmPtr->pick(brs)];
  auto lType = entry.decayChannels.at(prods).lType;

  if (lType == 0) {
    infoPtr->errorMsg("Warning in HadronWidths::pickDecay: "
      "channel does not set its angular momentum for ",
      std::to_string(idDec) + " -> " + std::to_string(prods.first) + " + " + std::to_string(prods.second));
    return false;
  }

  double mA, mB;
  if (!_pickMasses(prods.first, prods.second, m, mA, mB, lType))
    return false;

  idAOut = isAnti ? particleDataPtr->antiId(prods.first)  : prods.first;
  idBOut = isAnti ? particleDataPtr->antiId(prods.second) : prods.second;
  mAOut = mA; mBOut = mB;
  return true;
}


bool HadronWidths::pickExcitation(int idA, int idB, double eCM, 
  int& idCOut, double& mCOut, int& idDOut, double& mDOut) {

  // Pick an excitation channel
  vector<double> sigmas(excitationChannels.size());
  for (size_t i = 0; i < sigmas.size(); ++i)
    sigmas[i] = excitationChannels[i].sigma(eCM);
  auto& channel = excitationChannels[rndmPtr->pick(sigmas)];

  // The two nucleons have equal chance of becoming excited
  int maskA = channel.maskA, maskB = channel.maskB;  
  if (rndmPtr->flat() > 0.5)
    swap(maskA, maskB);

  // Construct ids of resonances from masks plus incoming ids
  int idCtmp = maskA + (idA - idA % 10);
  int idDtmp = maskB + (idB - idB % 10);
  
  // Pick masses
  double mCtmp, mDtmp;
  if (!_pickMasses(idCtmp, idDtmp, eCM, mCtmp, mDtmp, 1))
    return false;

  // Set output values and return
  idCOut = idCtmp; mCOut = mCtmp;
  idDOut = idDtmp; mDOut = mDtmp;
  return true;
}

}