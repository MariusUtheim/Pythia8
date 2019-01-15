
#include "Pythia8/LowEnergyData.h"

namespace Pythia8 {

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

bool LowEnergyData::readXML(istream& stream) {

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


double LowEnergyData::massDependentWidth(int id, double m) const {
  auto gen = speciesToGenus(id);
  return particleWidthPtr->width(gen, m);
}

double LowEnergyData::getStrangeness(int id) const {
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


double LowEnergyData::getStrangenessFactor(int id) const {
  int strangeness = getStrangeness(id);

  return strangeness / double(particleDataPtr->isBaryon(id) ? 3 : 2);
}


int LowEnergyData::getIsospin(int species) const {
  auto genus = speciesToGenus(species);
  return isospinType(genus);
}

int LowEnergyData::getIso3(int species) const {
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

double LowEnergyData::getClebschGordan2(int j1, int m1, int j2, int m2, int j, int m) const {
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

double LowEnergyData::getBR(int iRes, int i1, int i2, double eCM) const {
  ResGenus resR = speciesToGenus(iRes),
           resA = speciesToGenus(i1), resB = speciesToGenus(i2);
  
  double br = particleWidthPtr->branchingRatio(resR, resA + " " + resB, eCM);
  return br;
}


vector<int> LowEnergyData::getResonanceCandidates(int idA, int idB) const {
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

void LowEnergyData::print() const {

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

bool LowEnergyData::sanityCheck() {
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