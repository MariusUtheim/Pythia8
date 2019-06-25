
#include "Pythia8/ParticleWidths.h"

namespace Pythia8 {

typedef pair<int, int> keyType;

void ParticleWidths::Entry::addProducts(keyType key, Interpolator brs,
  int idA, int idB, int lType) {
    decayChannels.emplace(key, Channel(brs, idA, idB, lType));
}

double ParticleWidths::Entry::getWidth(keyType key, double eCM) const {
    auto iter = decayChannels.find(key);
    return (iter != decayChannels.end()) ? iter->second.br(eCM) * widths(eCM) : 0.;
}

double ParticleWidths::Entry::getBR(keyType key, double eCM) const {
    auto iter = decayChannels.find(key);
    return (iter != decayChannels.end()) ? iter->second.br(eCM) : 0.;
}

int ParticleWidths::Entry::getlType(keyType key) const {
    auto iter = decayChannels.find(key);
    return (iter != decayChannels.end()) ? iter->second.lType : 0;
}



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
keyType ParticleWidths::getKey(int& idR, int idA, int idB) const {

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


bool ParticleWidths::init(Info* infoPtrIn, Rndm* rndmPtrIn,
  ParticleData* particleDataPtrIn, string path) {
  infoPtr = infoPtrIn;
  rndmPtr = rndmPtrIn;
  particleDataPtr = particleDataPtrIn;

  ifstream stream(path);
  if (!stream.is_open()) {
    infoPtr->errorMsg( "Error in ParticleWidths::init: "
        "unable to open file");
    return false;
  }
  return readXML(stream);
}

bool ParticleWidths::readXML(istream& stream) {

  string line;

  while (getline(stream, line)) {

    string word1;
    istringstream(line) >> word1;

    if (word1 == "<width") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");
      if (id < 0) {
        infoPtr->errorMsg( "Error in ParticleWidths::readXML: "
          "Got negative id", std::to_string(id));
        continue;
      }
      auto* entry = particleDataPtr->findParticle(id);
      if (!entry) {
        infoPtr->errorMsg( "Error in ParticleWidths::readXML: "
          "Particle is not defined", std::to_string(id));
        continue;
      }
      if (!entry->isHadron()) {
        infoPtr->errorMsg( "Error in ParticleWidths::readXML: "
          "Particle is not defined as hadron", std::to_string(id));
        continue;
      }
      if (particleDataPtr->heaviestQuark(id) > 3) {
        infoPtr->errorMsg("Error in ParticleWidths:readXML: "
          "Particle contains a charmed or bottom hadron",
          std::to_string(id));
        continue;
      }

      // @TODO: Add a shift to avoid the boundary region on the left side
      double left = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");
      double m0 = doubleAttributeValue(line, "m0");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      this->entries.emplace(id, Entry(m0, Interpolator(left, right, data)));
    }
    else if (word1 == "<br") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");
      auto iter = entries.find(id);
      if (iter == entries.end()) {
        infoPtr->errorMsg( "Error in ParticleWidths::readXML: "
          "got br for a particle with undefined or ill-defined width");
        continue;
      }

      int lType = intAttributeValue(line, "lType");
      if (lType == 0) {
        infoPtr->errorMsg( "Error in ParticleWidths::readXML: "
          "lType is not defined");
        lType = 1;
      }

      istringstream productStr(attributeValue(line, "products"));
      vector<int> products; // @TODO this better
      bool gotInvalidParticle = false;
      int currentProduct;
      while (productStr >> currentProduct) {
        if (!particleDataPtr->isParticle(currentProduct)) gotInvalidParticle = true;
        products.push_back(currentProduct);
      }
      if (gotInvalidParticle) {
        infoPtr->errorMsg( "Error in ParticleWidths::readXML: "
          "decay product is not a particle",
          std::to_string(id) + " --> " + productStr.str());
        continue;
      }

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      keyType key = getKey(id, products[0], products[1]);
      
      auto& entry = iter->second;
      Interpolator br(entry.widths.left(), entry.widths.right(), data);
      entry.addProducts(key, br, products[0], products[1], lType);
    }
  }

  return true;
}



vector<int> ParticleWidths::getResonances() const {
  vector<int> resonances;
  for (auto& p : entries)
    resonances.push_back(p.first);
  return resonances;
}

double ParticleWidths::width(int id, double eCM) const {
  if (id < 0) id = -id;
  auto iter = entries.find(id);
  return (iter != entries.end()) ? iter->second.widths(eCM) : 0.;
}

double ParticleWidths::partialWidth(int idR, int idA, int idB, double m) const {
  auto prods = getKey(idR, idA, idB);
  auto iter = entries.find(idR);
  return (iter != entries.end()) ? iter->second.getWidth(prods, m) : 0.;
}

double ParticleWidths::branchingRatio(int idR, int idA, int idB, double m) const {
  auto prods = getKey(idR, idA, idB);
  auto iter = entries.find(idR);
  return (iter != entries.end()) ? iter->second.getBR(prods, m) : 0.;
}

double ParticleWidths::resonanceSigma(int idA, int idB, int idR,
  double eCM) const {
  
  // Ensure canonical ordering
  keyType key = getKey(idR, idA, idB);

  // Get width entry
  auto entryIter = entries.find(idR);
  if (entryIter == entries.end())
    return 0.;
  auto& entry = entryIter->second;

  // Find decay channel
  auto channelIter = entry.decayChannels.find(key);
  if (channelIter == entry.decayChannels.end())
    return 0.;
  auto& channel = channelIter->second;

  // Find particle entries
  auto* entryR = particleDataPtr->findParticle(idR);
  auto* entryA = particleDataPtr->findParticle(channel.idA);
  auto* entryB = particleDataPtr->findParticle(channel.idB);
  if (!entryR || !entryA || !entryB) {
    infoPtr->errorMsg("In ParticleWidths::resonanceSigma: "
      "got invalid particle id");
    return 0.;
  }

  // Get mass dependent width and branching ratios
  double gammaR = entry.widths(eCM);
  if (gammaR == 0) return 0.;

  double br = channel.br(eCM);
  if (br == 0) return 0.;

  // Calculate the resonance sigma
  double s = pow2(eCM), mA = entryA->m0(), mB = entryB->m0();
  double pCMS2 = 1 / (4 * s) * (s - pow2(mA + mB)) * (s - pow2(mA - mB));

  return GEVINVSQ2MB * M_PI / pCMS2
    * entryR->spinType() / (entryA->spinType() * entryB->spinType())
    * br * pow2(gammaR) / (pow2(entryR->m0() - eCM) + 0.25 * pow2(gammaR));
}


bool ParticleWidths::pickMasses(int idA, int idB, double eCM, double& mAOut, double& mBOut) {
  return _pickMasses(idA, idB, eCM, mAOut, mBOut, 1);
}

bool ParticleWidths::_pickMasses(int idA, int idB, double eCM,
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

static constexpr double MAX_LOOPS = 200;

bool ParticleWidths::_pickMass1(int idRes, double eCM, double mB, int lType,
  double& mAOut) {

  // Ensure resonance is positive - the mass distribution doesn't change
  idRes = abs(idRes);

  // Get width entry
  auto iter = entries.find(idRes);
  if (iter == entries.end()) {
    infoPtr->errorMsg("Error in ParticleWidths::pickMass: "
      "mass distribution for particle is not defined",
      std::to_string(idRes));
    return false;
  }
  Entry& entry(iter->second);

  // @TODO: Maybe an mPeak that is different from m0 will be more efficient
  double mMin = entry.widths.left(), 
         mMax = min(entry.widths.right(), eCM - mB),
         mPeak = particleDataPtr->m0(idRes),
         m0 = particleDataPtr->m0(idRes);

  // This can happen due to interpolation imprecision if eCM - mB is near mMin
  if (mMax < mMin)
    return false;

  if (mMin > m0 || mMin > mPeak) {
    infoPtr->errorMsg("Error in ParticleWidths::pickMass: "
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

  infoPtr->errorMsg("Warning in ParticleWidths::pickMass: "
    "Could not pick mass within prescribed number of iterations. ",
    std::to_string(idRes) + " in (" + std::to_string(mMin) + ", " + std::to_string(mMax) + ")");
  
  mAOut = mMin + rndmPtr->flat() * (mMax - mMin);
  return true;
}

// @TODO: Implement this properly
//        For the first iteration, we just pick id2 on shell
bool ParticleWidths::_pickMass2(int id1, int id2, double eCM, int lType,
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



bool ParticleWidths::pickDecay(int idDec, double m, int& idAOut, int& idBOut,
    double& mAOut, double& mBOut) {

  bool isAnti = (idDec < 0);
  if (isAnti)
    idDec = -idDec;

  auto entriesIter = entries.find(idDec);
  if (entriesIter == entries.end()) {
    infoPtr->errorMsg("Error in ParticleWidths::pickDecay: "
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
    infoPtr->errorMsg("Warning in ParticleWidths::pickDecay: "
      "no channels have positive branching ratios");
    return false;
  }
  
  auto prods = prodss[rndmPtr->pick(brs)];
  auto lType = entry.decayChannels.at(prods).lType;

  if (lType == 0) {
    infoPtr->errorMsg("Warning in ParticleWidths::pickDecay: "
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








}