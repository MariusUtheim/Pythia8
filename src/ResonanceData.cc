
#include "Pythia8/ResonanceData.h"


namespace Pythia8 {

double ResonanceData::getStrangeness(int id) const {
  // In baryon id xxxabcx, count number of occurrences of '3' in abc
  id = abs(id);
  int count = 0;
  for (int n = (id / 10) % 1000; n > 0; n /= 10) {
    int j = n % 10;
  
    if (j == 3)
      count++;
  }

  return count;
}

static string attributeValue(string line, string attribute) {

  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

/*
static bool boolAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  string tagLow = toLower(valString);
  return tagLow == "true" || tagLow == "1" || tagLow == "on"
      || tagLow == "yes" || tagLow == "ok";
}

static int intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;
}
*/

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

bool ResonanceData::readXML(istream& stream) {

  string line;

  while (getline(stream, line)) {

    string word1;
    istringstream(line) >> word1;

    if (word1 == "<resonanceClass") {
      completeTag(stream, line);

      ResClass classIdentifier = attributeValue(line, "identifier");

      istringstream generaStr(attributeValue(line, "genera"));
      vector<ResGenus> genera;
      ResGenus currentGenus;
      while (generaStr >> currentGenus) {
        genera.push_back(currentGenus);
        this->_genusToClass.emplace(currentGenus, classIdentifier);
      }
      
      this->_classToGenera.emplace(classIdentifier, genera);

    }
    else if (word1 == "<resonanceGenus") {
      completeTag(stream, line);

      ResGenus genusIdentifier = attributeValue(line, "identifier");

      istringstream speciesStr(attributeValue(line, "species"));
      vector<int> species;
      int currentSpecies;
      while (speciesStr >> currentSpecies) {
        species.push_back(currentSpecies);
        this->_speciesToGenus.emplace(currentSpecies, genusIdentifier);
      }

      this->_genusToSpecies.emplace(genusIdentifier, species);

      int strangeness = getStrangeness(species[0]);
      if (particleDataPtr->isMeson(species[0])) {
        if (strangeness == 0)
          this->_isospinType.emplace(genusIdentifier, 2 * (species.size() - 1));
        else if (strangeness == 1)
          this->_isospinType.emplace(genusIdentifier, 1);
        else if (strangeness == 2)
          this->_isospinType.emplace(genusIdentifier, 0);
        else
          throw " in ResonanceData::readXML: Strangeness does not give a sensible value"; // @TODO: Remove when we're confident getStrangeness works well
      }
      else {
        if (strangeness == 0)
          this->_isospinType.emplace(genusIdentifier, species.size() - 1);
        else if (strangeness == 1)
          this->_isospinType.emplace(genusIdentifier, species.size() == 1 ? 0 : 2);
        else if (strangeness == 2)
          this->_isospinType.emplace(genusIdentifier, 2);
        else if (strangeness == 3)
          this->_isospinType.emplace(genusIdentifier, 0);
        else
          throw " in ResonanceData::readXML: Strangeness does not give a sensible value"; // @TODO: Remove when we're confident getStrangeness works well
      }
    }
    else if (word1 == "<scatter") {
      completeTag(stream, line);

      istringstream inStr(attributeValue(line, "in"));
      pair<ResClass, ResClass> in;
      inStr >> in.first; inStr >> in.second;
      
      istringstream resonancesStr(attributeValue(line, "outs"));
      vector<pair<ResClass, ResClass>> outs;
      ResClass first;
      while (resonancesStr >> first) {
        pair<ResClass, ResClass> currentChannel;
        currentChannel.first = first;
        resonancesStr >> currentChannel.second;
        outs.push_back(currentChannel);
      }
      
      this->scatterChannels.emplace(in, outs);

    }
    else if (word1 == "<sigmaTotal") {
      completeTag(stream, line);

      istringstream inStr(attributeValue(line, "in"));
      pair<ResClass, ResClass> in;
      inStr >> in.first; inStr >> in.second;

      double left = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      this->totalSigmaDistribution.emplace(in, Interpolator(left, right, data));
    }
    else if (word1 == "<brFactor") {
      completeTag(stream, line);

      istringstream outStr(attributeValue(line, "out"));
      pair<ResGenus, ResGenus> out;
      outStr >> out.first; outStr >> out.second;

      double left = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      this->partialSigmaDistribution.emplace(out, Interpolator(left, right, data));
    }
  }

  return true;
}

//--------------------------------------------------------------------------

double ResonanceData::getTotalSigma(int idA, int idB, double eCM) const {

  if (eCM < particleDataPtr->m0(idA) + particleDataPtr->m0(idB))
    return 0.;

  if (particleDataPtr->isBaryon(idA) && particleDataPtr->isBaryon(idB)) {
    if (idA * idB < 0) {
      if (idA > 0)
        return getTotalSigmaBBbar(idA, idB, eCM);
      else
        return getTotalSigmaBBbar(idB, idA, eCM);
    }
    else
      return getTotalSigmaBB(idA, idB, eCM);
  }
  else {
    if (particleDataPtr->isMeson(idB))
      return getTotalSigmaXM(idA, idB, eCM);
    else
      return getTotalSigmaXM(idB, idA, eCM);
  }
  
  return 0.;
}


//--------------------------------------------------------------------------

// B+B section

double ResonanceData::getTotalSigmaBB(int idA, int idB, double eCM) const {
  cout << "ResonanceData::getTotalSigmaBB not implemented" << endl;
  return 0.;
}

const Interpolator& ResonanceData::getDiffractiveSigmaDistribution(int idA, int idB) const {
  auto cls = classify(idA, idB);
  auto ptr = totalSigmaDistribution.find(cls);
  if (ptr == totalSigmaDistribution.end())
    return Interpolator::Zero;
  else
    return ptr->second;
}

double ResonanceData::getDiffractiveSigma(int idA, int idB, double eCM) const {
  auto cls = classify(idA, idB);
  auto ptr = totalSigmaDistribution.find(cls);
  if (ptr == totalSigmaDistribution.end())
    return 0.;
  else
    return ptr->second(eCM);
}

vector<pair<int, int>> ResonanceData::getDiffractiveOutputs(int idA, int idB) const {
  ResClass clsA = speciesToClass(idA), clsB = speciesToClass(idB);
  auto cls = pair<ResClass, ResClass>(clsA, clsB);

  int totalCharge = particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB);
  auto& outputClasses = scatterChannels.at(cls);

  vector<pair<int, int>> possibleOutputs;

  for (auto outputClass : outputClasses) {

    auto parCs = classToSpecies(outputClass.first), 
         parDs = classToSpecies(outputClass.second);

    if (idA > 0 && idB > 0) {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge) {
            possibleOutputs.push_back(pair<int, int>(pC, pD));

          }
    }
    else if (idA < 0 && idB < 0) {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge)
            possibleOutputs.push_back(pair<int, int>(-pC, -pD));
    }
    else {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge) {
            possibleOutputs.push_back(pair<int, int>(-pC, pD));
            possibleOutputs.push_back(pair<int, int>(pC, -pD));
          }
    }

  }
  
  return possibleOutputs;
}

vector<pair<pair<int, int>, double>> ResonanceData::getOutputsWithFrequencies(int idA, int idB, double eCM) const {
  
  auto outs = getDiffractiveOutputs(idA, idB);
  vector<pair<pair<int, int>, double>> outsWithFreq(outs.size());

  for (int iR = 0; iR < (int)outs.size(); ++iR) {
    auto gens = genify(outs[iR].first, outs[iR].second);
    double weight = partialSigmaDistribution.at(gens)(eCM);
    outsWithFreq[iR] = pair<pair<int, int>, double>(outs[iR], weight);
  }

  return outsWithFreq;
}

//--------------------------------------------------------------------------

// B+Bbar section

double ResonanceData::getTotalSigmaBBbar(int idA, int idB, double eCM) const {
  cout << "ResonanceData::getTotalSigmaBBbar not implemented" << endl;
  return 0.;
}

double ResonanceData::getStrangenessFactor(int id) const {
  int strangeness = getStrangeness(id);

  return strangeness / double(particleDataPtr->isBaryon(id) ? 3 : 2);
}

// Calculate sigma_annihilation. Assumes one index refers to a baryon, the other to an antibaryon
double ResonanceData::getAnnihilationSigma(int idA, int idB, double eCM) const {
  if (idA > idB) 
    swap(idA, idB);

  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM <= mA + mB)
    return 0.;

  double xsA = getStrangenessFactor(idA), xsB = getStrangenessFactor(idB);

  static const double sigma0N = 120., s0 = 3.5214, A2 = 0.0025, B = 0.6;
  double s = eCM * eCM;
  double sigmaNN = sigma0N * s0/s * (A2 * s0 / (pow2(s - s0) + A2 * s0) + B);

  // @TODO separate into annihilation, diffractive and elastic

  return sigmaNN * (1. - 0.4 * xsA) * (1. - 0.4 * xsB);
}

//--------------------------------------------------------------------------

// X+M section

static const vector<double> ppiminusInterpolData = 
  { 20.3919, 35.1056, 53.4896, 67.8302, 61.3948, 47.1339, 36.1912, 
    29.9881, 26.8335, 26.5702, 27.7574, 28.5599, 30.0137, 32.8038, 
    37.3226, 44.0601, 45.7762, 42.1849, 37.6042, 35.9693, 38.7317, 
    45.5506, 53.6462, 58.0236, 54.3861, 47.0179, 41.6299, 38.2877, 
    36.5429, 36.1408, 36.1086, 36.3178, 36.1149, 35.7958, 35.3629, 
    34.9317, 34.5016, 34.3947, 34.5755, 34.7563, 34.9371, 35.1311, 
    35.3531, 35.575, 35.797, 35.9509, 35.7146, 35.4783, 35.2419, 35.0056, 
    34.7693, 34.5329, 34.2954, 34.0514, 33.8073, 33.5632, 33.3192, 
    33.0751, 32.831, 32.587, 32.381, 32.278, 32.1749, 32.0718, 31.9688, 
    31.8657, 31.7626, 31.6595, 31.5565, 31.4534, 31.3503, 31.244, 
    31.1329, 31.0218, 30.9108, 30.7997, 30.6887, 30.5776, 30.4666, 
    30.3555, 30.2445, 30.1334, 30.0257, 29.9217, 29.8176, 29.7135, 
    29.6095, 29.5054, 29.4014, 29.2973, 29.1932, 29.0978, 29.0211, 
    28.9445, 28.8678, 28.7911, 28.7144, 28.6378, 28.5611, 28.4844, 
    28.4077, 28.3311, 28.2806, 28.2432, 28.2058, 28.1684, 28.131, 
    28.0936, 28.0562, 28.0188, 27.9814, 27.944, 27.9066, 27.8692, 
    27.8318, 27.7944, 27.757, 27.7196, 27.6822, 27.6448, 27.6074 };

static const Interpolator ppiminusInterpol(1.15, 3.89, ppiminusInterpolData);

static const vector<double> ppiplusInterpolData =
 { 33.1408, 30.3605, 29.3346, 29.8405, 30.06, 30.1206, 29.917, 29.3581,
   28.8239, 28.5466, 28.2694, 27.9941, 27.7218, 27.4494, 27.1523,
   26.8496, 26.547, 26.3317, 26.1572, 25.9828, 25.8343 };
static const Interpolator ppiplusInterpol(2.03, 3.50, ppiplusInterpolData);

double ResonanceData::getTotalSigmaXM(int idX, int idM, double eCM) const {
  if (idX == 2212 && idM == -211) {
    return eCM <= 1.8 ? getResonanceSigma(idX, idM, eCM)
         : eCM <= 3.5 ? ppiminusInterpol(eCM)
         :              13.7 * pow(eCM, 0.158) + 35.9 * pow(eCM, -0.90);
  }
  else if (idX == 2212 && idM == 211) {
    return eCM <= 2.03 ? getResonanceSigma(idX, idM, eCM)
         : eCM <= 3.50 ? ppiplusInterpol(eCM)
         :              13.7 * pow(eCM, 0.158) + 35.9 * pow(eCM, -0.90);
  }
  else
    return 0.;
}

void ResonanceData::showPickProbabilities(int idX, int idM, double eCM) const {
  if (idX == 2212 && idM == -211) {
    double total = getTotalSigmaXM(idX, idM, eCM);
    
    double resonant = getResonanceSigma(idX, idM, eCM);
    
    double mp = particleDataPtr->m0(2212), mpi = particleDataPtr->m0(211);
    double pLab = sqrt((pow2(eCM) - pow2(mp + mpi)) * (pow2(eCM) - (mp - mpi)) / (2 * mp));
    double elastic = 1.76 + 11.2 * pow(pLab, -0.64) + 0.043 * pow2(log(pLab))
                   - getElasticResonanceSigma(idX, idM, eCM);

    double string = total - resonant - elastic;

    cout << setw(8) << "Total: " << total << endl
         << setw(8) << "Res: " << resonant << endl
         << setw(8) << "El: " << elastic << endl
         << setw(8) << "Str: " << string << endl
         << endl;
  
  }
}

double ResonanceData::getBR(int idR, int idA, int idB, double eCM) const {
  ResGenus resR = speciesToGenus(idR),
    resA = speciesToGenus(idA), resB = speciesToGenus(idB);
  
  double br = particleWidthPtr->branchingRatio(resR, resA + " " + resB, eCM);
  return br;
}

int ResonanceData::getIsospin(int species) const {
  auto genus = speciesToGenus(species);
  return isospinType(genus);
}

int ResonanceData::getIso3(int species) const {
  int totalSpin = getIsospin(species);
  int charge = particleDataPtr->chargeType(abs(species)) / 3;
  int sign = species > 0 ? 1 : -1;
  return sign * 
       ( totalSpin == 0 ? 0
       : totalSpin == 1 ? 2 * charge - 1
       : totalSpin == 2 ? 2 * charge 
       : totalSpin == 3 ? 2 * charge - 1
       : 0 );
}

double ResonanceData::getClebschGordan2(int j1, int m1, int j2, int m2, int j, int m) const {
  // @TODO: Maybe don't give such aggressive warning messages. For now these are kept in as
  //        tests, as in the current code, they should never be hit, but maybe allow it in the future
  if (j1 < 0 || j2 < 0 || j < 0)
  { cout << "ResonanceData::getClebschGordan: Got negative isospins" << endl; return NAN; }
  if (j1 > 3 || j2 > 3)
  { cout << "ResonanceData::getClebschGordan: Got isospin greater than 3/2; case is unhandled" << endl; return NAN; }
  if (m1 < -j1 || m1 > j1 || m2 < -j2 || m2 > j2 || m < -j || m > j)
  { cout << "ResonanceData::getClebschGordan: component 3 out of range" << endl; return NAN; }
  if ((j1 + j2 - j) % 2 != 0)
  { cout << "ResonanceData::getClebschGordan: integer/half-integerness is not conserved" << endl; return NAN; }
  if (m != m1 + m2)
  { cout << "ResonanceData::getClebschGordan: component 3 not conserved" << endl; return NAN; }
  if (j < abs(j1 - j2) || j > j1 + j2)
  { cout << "ResonanceData::getClebschGordan: outgoing isospin out of range" << endl; return NAN; }
  if ((j1 + m1) % 2 != 0 || (j2 + m2) % 2 != 0 || (j + m) % 2 != 0)
  { cout << "ResonanceData::getClebschGordan: integer-spin particle has half-spin component 3, or vice versa" << endl; return NAN; }

  // Cover all boundary cases
  if (j1 + j2 == j && abs(m1 + m2) == j)
    return 1.;

  if (j1 > j2)
  { swap(j1, j2); swap(m1, m2); }

  // We have already checked that spin values are within physical range,
  // so in this case it is implicit that m1 = 0, j2 = j and m2 = m.
  if (j1 == 0)
    return 1.;

  // We have checked all boundary cases, so here m1 + m2 = 0
  if (j1 == 1 && j2 == 1)
    return 1. / 2.;
  else if (j1 == 1 && j2 == 2) {
    if (abs(m2 + m) == j)
      return 1. / 3.;
    else
      return 2. / 3.;
  }
  else if (j1 == 1 && j2 == 3) {
    if (m == 0)
      return 1. / 2.;
    else if ((j == 2 && m1 * m2 > 0) || (j == 4 && m1 * m2 < 0))
      return 1. / 4.;
    else
      return 3. / 4.;
  }
  else if (j1 == 2 && j2 == 2) {
    if (m != 0)
      return 1. / 2.;
    else if (m1 == 0 && m2 == 0) {
      if (j == 0) return 1. / 3.;
      else if (j == 2) return 0.;
      else return 2. / 3.;
    }
    else {
      if (j == 0) return 1. / 3.;
      else if (j == 2) return 1. / 2.;
      else return 1. / 6.;
    }
  }
  else if (j1 == 2 && j2 == 3) {
    if (abs(m) == 3) 
      return (j + abs(m1) == 5) ? 2./5. : 3./5.;
    else {
      static const double table[3][3] = {
        { 1./2., 2./5., 1./10. },
        { 1./3., 1./15., 3./5. },
        { 1./6., 8./15., 3./10. }
      };
      int mtype = 1 + (m > 0 ? m1/2 : -m1/2);
      int jtype = (j - 1) / 2;
      return table[mtype][jtype];
    }
  }
  else if (j1 == 3 && j2 == 3) {
    // @TODO: Implement this
    cout << "ResonanceData::getClebschGordan: 3/2 x 3/2 not implemented" << endl;
    return 1.;
  }
  else {
    cout << "ResonanceData::getClebschGordan: got spins higher than 3/2" << endl;
    return 1.;
  }
}

vector<int> ResonanceData::getResonanceCandidates(int idA, int idB) const {
  auto gens = genify(abs(idA), abs(idB));

  auto totalCharge = particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB);
  auto totalIso3 = getIso3(idA) + getIso3(idB);

  vector<int> resonanceCandidates;

  if (particleDataPtr->isMeson(idA)) {
    for (auto genus : classToGenera("M*")) {
      for (auto species :  genusToSpecies(genus))
        if (particleDataPtr->chargeType(species) == totalCharge
        &&  getIso3(species) == totalIso3) 
          resonanceCandidates.push_back(species); // break; 
    }
  }
  else {
    int strangeness = abs(getStrangeness(idA) + getStrangeness(idB));
    auto classes = strangeness == 0 ? vector<string>({ "N*", "D*", "D" })
                 : strangeness == 1 ? vector<string>({ "S1" })
                 : strangeness == 2 ? vector<string>({ "S2" })
                 : vector<string>();
    for (auto cls : classes)
    for (auto genus : classToGenera(cls)) {
      for (auto species : genusToSpecies(genus)) {
        if (idA < 0) 
          species = -species;
        if (particleDataPtr->chargeType(species) == totalCharge
        &&  getIso3(species) == totalIso3)
          resonanceCandidates.push_back(species); // break; 
      }
    }
  }

  return resonanceCandidates;
}

double ResonanceData::getPartialResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const {
  double br = getBR(idR, idA, idB, eCM);
  if (br == 0.)
    return 0.;

  double gammaRes2 = pow2(particleWidthPtr->mass(speciesToGenus(idR), eCM));

  int iA = getIsospin(idA);
  int i3A = getIso3(idA);

  int iB = getIsospin(idB);
  int i3B = getIso3(idB);
  
  int iR = getIsospin(idR);
  int i3R = getIso3(idR);

  double cg2 = getClebschGordan2(iA, i3A, iB, i3B, iR, i3R);

  if (isnan(cg2)) {
    cout << " for " << particleDataPtr->name(idA) << " + " << particleDataPtr->name(idB) << " --> " << particleDataPtr->name(idR) << endl
         << "     < " << iA << " " << i3A << " , " << iB << " " << i3B << " | " << iR << " " << i3R << " > " << endl;
  }

  int nJRes = particleDataPtr->spinType(idR);
  double m0 = particleDataPtr->m0(idR);

  double contribution = cg2 * nJRes * br * gammaRes2/(pow2(eCM - m0) + 0.25 * gammaRes2);

  if (gensEqual && idA != idB)
    contribution *= 2;

  return contribution;
}

double ResonanceData::getPartialElasticResonanceSigma(int idA, int idB, int idR, bool gensEqual, double eCM) const {
  double br = getBR(idR, idA, idB, eCM);
  if (br == 0.)
    return 0.;

  double gammaRes2 = pow2(particleWidthPtr->mass(speciesToGenus(idR), eCM));

  int iA = getIsospin(idA);
  int i3A = getIso3(idA);

  int iB = getIsospin(idB);
  int i3B = getIso3(idB);
  
  int iR = getIsospin(idR);
  int i3R = getIso3(idR);

  double cg2 = getClebschGordan2(iA, i3A, iB, i3B, iR, i3R);

  if (isnan(cg2)) {
    cout << " for " << particleDataPtr->name(idA) << " + " << particleDataPtr->name(idB) << " --> " << particleDataPtr->name(idR) << endl
         << "     < " << iA << " " << i3A << " , " << iB << " " << i3B << " | " << iR << " " << i3R << " > " << endl;
  }

  int nJRes = particleDataPtr->spinType(idR);
  double m0 = particleDataPtr->m0(idR);

  double contribution = cg2 * nJRes * br * gammaRes2/(pow2(eCM - m0) + 0.25 * gammaRes2);

  if (gensEqual && idA != idB)
    contribution *= 2;

  contribution *= br * cg2; // Elastic correction

  return contribution;
}

// @TODO Probably take Particles instead of just ids
double ResonanceData::getResonanceSigma(int idA, int idB, double eCM) const {

  // For K_S and K_L, take average of K and Kbar
  if (idA == 310 || idA == 130)
    return 0.5 * (getResonanceSigma(311, idB, eCM) + getResonanceSigma(-311, idB, eCM));
  if (idB == 310 || idB == 130)
    return 0.5 * (getResonanceSigma(idA, 311, eCM) + getResonanceSigma(idA, -311, eCM));

  // @TODO Deal with what happens when there is an antibaryon
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM < mA + mB) return 0.; // @TODO: This isn't needed if we take two input particles and calculate s

  // Ensure baryon is always idA, if there is a baryon
  if (particleDataPtr->isBaryon(idB) && particleDataPtr->isMeson(idA))
    swap(idA, idB);


  vector<int> resonanceCandidates = getResonanceCandidates(idA, idB);
  bool gensEqual = speciesToGenus(abs(idA)) == speciesToGenus(abs(idB));

  double sigmaRes = 0;
  for (auto idR : resonanceCandidates) {
    sigmaRes += getPartialResonanceSigma(idA, idB, idR, gensEqual, eCM);
  }

  double s = eCM * eCM;
  double pCMS2 = (s - pow2(mA + mB)) * (s - pow2(mA - mB)) / (4 * s);

  // @TODO define constant 0.38937966 = GeV^-2 to mb
  double preFactor = 0.38937966 * M_PI / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB) * pCMS2);
  sigmaRes *= preFactor;

  return sigmaRes;
}

double ResonanceData::getElasticResonanceSigma(int idA, int idB, double eCM) const {

  // For K_S and K_L, take average of K and Kbar
  if (idA == 310 || idA == 130)
    return 0.5 * (getResonanceSigma(311, idB, eCM) + getResonanceSigma(-311, idB, eCM));
  if (idB == 310 || idB == 130)
    return 0.5 * (getResonanceSigma(idA, 311, eCM) + getResonanceSigma(idA, -311, eCM));

  // @TODO Deal with what happens when there is an antibaryon
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM < mA + mB) return 0.; // @TODO: This isn't needed if we take two input particles and calculate s

  // Ensure baryon is always idA, if there is a baryon
  if (particleDataPtr->isBaryon(idB) && particleDataPtr->isMeson(idA))
    swap(idA, idB);


  vector<int> resonanceCandidates = getResonanceCandidates(idA, idB);
  bool gensEqual = speciesToGenus(abs(idA)) == speciesToGenus(abs(idB));

  double sigmaRes = 0;
  for (auto idR : resonanceCandidates) {
    sigmaRes += getPartialElasticResonanceSigma(idA, idB, idR, gensEqual, eCM);
  }

  double s = eCM * eCM;
  double pCMS2 = (s - pow2(mA + mB)) * (s - pow2(mA - mB)) / (4 * s);

  // @TODO define constant 0.38937966 = GeV^-2 to mb
  double preFactor = 0.38937966 * M_PI / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB) * pCMS2);
  sigmaRes *= preFactor;

  return sigmaRes;
}



//--------------------------------------------------------------------------


void ResonanceData::print() const {

  cout << "== Classes ==" << endl;

  for (auto cls : _classToGenera) {
    cout << cls.first << " = { ";
    for (auto gen : cls.second)
      cout << gen << " ";
    cout << "}" << endl;
  }

  cout << endl << "== Genera ==" << endl;
  for (auto gen : _genusToSpecies) {
    cout << gen.first << "(" << getIsospin(gen.second[0]) << ") = { ";
    for (auto spc : gen.second)
      cout << spc << " ";
    cout << "}" << endl;
  }

  cout << endl << "== Particles in each class ==" << endl;
  for (auto cls : _classToGenera) {
    cout << cls.first << " = { ";

    for (auto gen : cls.second)
      for (auto spc : _genusToSpecies.at(gen))
        cout << spc << " ";
    cout << "}" << endl;
  }

  cout << endl;
}

bool ResonanceData::sanityCheck() {
  bool everythingWorks = true;

  for (auto cls : _classToGenera) {
    for (auto gen : cls.second) {
      if (_genusToSpecies.find(gen) == _genusToSpecies.end()) {
        cout << " ERROR: genus " << gen << " of class " << cls.first << " is missing" << endl;
        everythingWorks = false;
      }
    }
  }

  for (auto gen : _genusToSpecies) {
    for (auto spc : gen.second) {
      if (!particleDataPtr->isParticle(spc)) {
        cout << " ERROR: particle species " << spc << " of genus " << gen.first << " does not exist" << endl;
        everythingWorks = false;
      }
    }

    if (_genusToClass.find(gen.first) == _genusToClass.end()) {
      cout << " WARNING: genus " << gen.first << " does not belong to any class" << endl;
    }
  }

  for (auto clsA : _classToGenera) 
  for (auto clsB : _classToGenera) {
    pair<ResClass, ResClass> clsPair(clsA.first, clsB.first);
    if (totalSigmaDistribution.find(clsPair) != totalSigmaDistribution.end()) {
      if (scatterChannels.find(clsPair) == scatterChannels.end()) {
        cout << " ERROR: Total sigma for " << clsPair.first << " " << clsPair.second << " is nonzero, but scattering channel not found" << endl;
        everythingWorks = false;
      }
      else
      {
        for (auto chn : scatterChannels.at(clsPair))
        for (auto genC : _classToGenera.at(chn.first))
        for (auto genD : _classToGenera.at(chn.second)) {
          pair<ResGenus, ResGenus> genPair(genC, genD);
          if (partialSigmaDistribution.find(genPair) == partialSigmaDistribution.end()) {
            cout << " ERROR: Total sigma for " << clsPair.first << " " << clsPair.second << " is nonzero, but relative frequency for output " << genPair.first << " " << genPair.second << " not found" << endl;
            everythingWorks = false;
          }
        }
      }
    }
  }
    
  return everythingWorks;
}

}
