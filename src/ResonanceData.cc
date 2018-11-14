
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
        this->genusToClass.emplace(currentGenus, classIdentifier);
      }
      
      this->classToGenera.emplace(classIdentifier, genera);

    }
    else if (word1 == "<resonanceGenus") {
      completeTag(stream, line);

      ResGenus genusIdentifier = attributeValue(line, "identifier");

      istringstream speciesStr(attributeValue(line, "species"));
      vector<int> species;
      int currentSpecies;
      while (speciesStr >> currentSpecies) {
        species.push_back(currentSpecies);
        this->speciesToGenus.emplace(currentSpecies, genusIdentifier);
      }
      
      this->genusToParticles.emplace(genusIdentifier, species);

      int strangeness = getStrangeness(species[0]);
      if (particleDataPtr->isMeson(species[0])) {
        if (strangeness == 0)
          this->isospinType.emplace(genusIdentifier, 2 * (species.size() - 1));
        else if (strangeness == 1)
          this->isospinType.emplace(genusIdentifier, 1);
        else if (strangeness == 2)
          this->isospinType.emplace(genusIdentifier, 0);
        else
          throw " in ResonanceData::readXML: Strangeness does not give a sensible value"; // @TODO: Remove when we're confident getStrangeness works well
      }
      else {
        if (strangeness == 0)
          this->isospinType.emplace(genusIdentifier, species.size() - 1);
        else if (strangeness == 1)
          this->isospinType.emplace(genusIdentifier, species.size() == 1 ? 0 : 2);
        else if (strangeness == 2)
          this->isospinType.emplace(genusIdentifier, 2);
        else if (strangeness == 3)
          this->isospinType.emplace(genusIdentifier, 0);
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

  double sigmaTotal = 0;

  if (particleDataPtr->isBaryon(idA) && particleDataPtr->isBaryon(idB)) {
    if (idA * idB < 0) 
      sigmaTotal += getAnnihilationSigma(idA, idB, eCM);
    else
      sigmaTotal += getDiffractiveSigma(idA, idB, eCM);
  }
  else
    sigmaTotal += getResonanceSigma(idA, idB, eCM);

  sigmaTotal += getElasticSigma(idA, idB, eCM);

  return sigmaTotal;
}


//--------------------------------------------------------------------------

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
  ResClass clsA = classify(idA), clsB = classify(idB);
  auto cls = classify(idA, idB);

  int totalCharge = particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB);
  auto& outputClasses = scatterChannels.at(cls);

  vector<pair<int, int>> possibleOutputs;

  for (auto outputClass : outputClasses) {

    //cout << clsA << " " << clsB << " --> " << outputClass.first << " " << outputClass.second << endl;

    auto parCs = classToSpecies(outputClass.first), 
         parDs = classToSpecies(outputClass.second);

    if (idA > 0 && idB > 0) {
      for (auto pC : parCs)
        for (auto pD : parDs)
          if (particleDataPtr->chargeType(pC) + particleDataPtr->chargeType(pD) == totalCharge) {
            possibleOutputs.push_back(pair<int, int>(pC, pD));

            //cout << setw(9) << left << particleDataPtr->name(pC) << setw(14) << particleDataPtr->name(pD);
            //auto gens = genify(pC, pD);
            //double weight = partialSigmaDistribution.at(gens)(5.);
            //cout << weight << endl;
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

  cout << " --- Annihilation between " << particleDataPtr->name(idA) << " and " << particleDataPtr->name(idB) << " --- " << endl
       << "sigmaNN = " << sigmaNN << endl
       << "xsA = " << xsA << endl
       << "xsB = " << xsB << endl
       << endl;
  

  return sigmaNN * (1. - 0.4 * xsA) * (1. - 0.4 * xsB);
}

double ResonanceData::getElasticSigma(int idA, int idB, double eCM) const {
  return 0. * idA * idB * eCM; // @TODO something
}

//--------------------------------------------------------------------------


double ResonanceData::getBR(int idR, int idA, int idB, double eCM) const {
  ResGenus resR = speciesToGenus.at(abs(idR)),
    resA = speciesToGenus.at(abs(idA)), resB = speciesToGenus.at(abs(idB));
  
  double br = particleWidthPtr->branchingRatio(resR, resA + " " + resB, eCM);
  return br;
}

int ResonanceData::getIsospin(int species) const {
  auto genus = speciesToGenus.at(abs(species));
  return isospinType.at(genus);
}

int ResonanceData::getIso3(int species) const {
  int totalSpin = getIsospin(species);
  int charge = particleDataPtr->chargeType(abs(species)) / 2;
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
  { cout << "ResonanceData::getClebschGordan: Got negative isospins" << endl; return 0.; }
  if (j1 > 3 || j2 > 3)
  { cout << "ResonanceData::getClebschGordan: Got isospin greater than 3/2; case is unhandled" << endl; return 0.; }
  if (m1 < -j1 || m1 > j1 || m2 < -j2 || m2 > j2 || m < -j || m > j)
  { cout << "ResonanceData::getClebschGordan: component 3 out of range" << endl; return 0.; }
  if ((j1 + j2 - j) % 2 != 0)
  { cout << "ResonanceData::getClebschGordan: parity is not conserved" << endl; return 0.; }
  if (m != m1 + m2)
  { cout << "ResonanceData::getClebschGordan: component 3 not conserved" << endl; return 0.; }
  if (j < abs(j1 - j2) || j > j1 + j2)
  { cout << "ResonanceData::getClebschGordan: outgoing isospin out of range" << endl; return 0.; }
  if ((j1 + m1) % 2 != 0 || (j2 + m2) % 2 != 0 || (j + m) % 2 != 0)
  { cout << "ResonanceData::getClebschGordan: integer-spin particle has half-spin component 3, or vice versa" << endl; return 0.; }

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

// @TODO Probably take Particles instead of just ids
double ResonanceData::getResonanceSigma(int idA, int idB, double eCM) const {
  // @TODO Deal with what happens when there is an antibaryon
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM < mA + mB) return 0.; // @TODO: This isn't needed if we take two input particles and calculate s

  // Ensure baryon is always idA, if there is a baryon
  if (particleDataPtr->isBaryon(idB) && particleDataPtr->isMeson(idA))
    swap(idA, idB);
  auto gens = genify(abs(idA), abs(idB));

  auto totalCharge = particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB);
  vector<int> resonanceCandidates;

  double sigmaRes = 0;

  if (particleDataPtr->isMeson(idA)) {
    for (auto genus : classToGenera.at("M*")) {
      for (auto species :  genusToParticles.at(genus))
        if (particleDataPtr->chargeType(species) == totalCharge) 
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
    for (auto genus : classToGenera.at(cls)) {
      for (auto species :  genusToParticles.at(genus))
        if (particleDataPtr->chargeType(species) == totalCharge) 
          resonanceCandidates.push_back(species); // break; 
    }
  }


  cout << " === Colliding " << particleDataPtr->name(idA) << " with " << particleDataPtr->name(idB) << " @ " << eCM << " GeV === " << endl << endl;

  for (auto idR : resonanceCandidates) {
    double br = getBR(idR, idA, idB, eCM);
    if (br == 0.) {
      cout << "Zero BR for " << speciesToGenus.at(abs(idR)) << endl << endl;
      continue;
    }

    double gammaRes2 = pow2(particleWidthPtr->mass(speciesToGenus.at(idR), eCM));


    int iA = getIsospin(idA);
    int i3A = getIso3(idA);

    int iB = getIsospin(idB);
    int i3B = getIso3(idB);
    
    int iR = getIsospin(idR);
    int i3R = getIso3(idR);
    
    double cg2 = getClebschGordan2(iA, i3A, iB, i3B, iR, i3R);


    int nJRes = particleDataPtr->spinType(idR);
    double m0 = particleDataPtr->m0(idR);

    double contribution = cg2 * nJRes * br * gammaRes2/(pow2(eCM - m0) + 0.25 * gammaRes2);

    if (gens.first == gens.second && idA != idB)
      contribution *= 2;

    cout << "Contribution from " << particleDataPtr->name(idR) << ": " << contribution << endl
          << " Clebsch-Gordan = <" << iA << " " << i3A << " , " << iB << " " << i3B << " | " << iR << " " << i3R << " >^2 = " << cg2 << endl
          << " (2S_R + 1) = " << nJRes << endl
          << " Branching ratio " << br << endl
          << " Gamma^2 = " << gammaRes2 << endl
          << " E_CM - m0 = " << eCM - m0 << endl
          << endl;

    sigmaRes += contribution;
  }

  double s = eCM * eCM;
  double pCMS2 = (s - pow2(mA + mB)) * (s - pow2(mA - mB)) / (4 * s);

  // @TODO define constant 0.38937966 = GeV^-2 to mb
  double preFactor = 0.38937966 * M_PI / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB) * pCMS2);
  cout << "Prefactor = " << preFactor << endl;
  sigmaRes *= preFactor;

  if (sigmaRes > 0 && particleDataPtr->isMeson(idA))
    sigmaRes += 5.;

  return sigmaRes;
}

//--------------------------------------------------------------------------


void ResonanceData::print() const {

  cout << "== Classes ==" << endl;

  for (auto cls : classToGenera) {
    cout << cls.first << " = { ";
    for (auto gen : cls.second)
      cout << gen << " ";
    cout << "}" << endl;
  }

  cout << endl << "== Genera ==" << endl;
  for (auto gen : genusToParticles) {
    cout << gen.first << "(" << getIsospin(gen.second[0]) << ") = { ";
    for (auto spc : gen.second)
      cout << spc << " ";
    cout << "}" << endl;
  }

  cout << endl << "== Particles in each class ==" << endl;
  for (auto cls : classToGenera) {
    cout << cls.first << " = { ";

    for (auto gen : cls.second)
      for (auto spc : genusToParticles.at(gen))
        cout << spc << " ";
    cout << "}" << endl;
  }

  cout << endl;
}

bool ResonanceData::sanityCheck() {
  bool everythingWorks = true;

  for (auto cls : classToGenera) {
    for (auto gen : cls.second) {
      if (genusToParticles.find(gen) == genusToParticles.end()) {
        cout << " ERROR: genus " << gen << " of class " << cls.first << " is missing" << endl;
        everythingWorks = false;
      }
    }
  }

  for (auto gen : genusToParticles) {
    for (auto spc : gen.second) {
      if (!particleDataPtr->isParticle(spc)) {
        cout << " ERROR: particle species " << spc << " of genus " << gen.first << " does not exist" << endl;
        everythingWorks = false;
      }
    }

    if (genusToClass.find(gen.first) == genusToClass.end()) {
      cout << " WARNING: genus " << gen.first << " does not belong to any class" << endl;
    }
  }

  for (auto clsA : classToGenera) 
  for (auto clsB : classToGenera) {
    pair<ResClass, ResClass> clsPair(clsA.first, clsB.first);
    if (totalSigmaDistribution.find(clsPair) != totalSigmaDistribution.end()) {
      if (scatterChannels.find(clsPair) == scatterChannels.end()) {
        cout << " ERROR: Total sigma for " << clsPair.first << " " << clsPair.second << " is nonzero, but scattering channel not found" << endl;
        everythingWorks = false;
      }
      else
      {
        for (auto chn : scatterChannels.at(clsPair))
        for (auto genC : classToGenera.at(chn.first))
        for (auto genD : classToGenera.at(chn.second)) {
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
