
#include "Pythia8/ResonanceData.h"


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
      if (particleDataPtr->isMeson(species[0]))
        this->isospinType.emplace(genusIdentifier, 2 * (species.size() - 1));
      else
        this->isospinType.emplace(genusIdentifier, species.size() - 1);
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

vector<pair<int, int>> ResonanceData::getPossibleOutputs(int idA, int idB) const {
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

//--------------------------------------------------------------------------


double ResonanceData::getBR(int idR, int idA, int idB) const {
  auto* dat = particleDataPtr->findParticle(idR); // @TODO: Validate output
  
  for (int iChannel = 0; iChannel < dat->sizeChannels(); ++iChannel) {
    auto channel = dat->channel(iChannel);
    if ((channel.product(0) == idA && channel.product(1) == idB)
    ||  (channel.product(0) == idB && channel.product(1) == idA))
      return channel.bRatio();
  }

  return 0.;
}

int ResonanceData::getIsospin(int species) const {
  auto genus = speciesToGenus.at(abs(species));
  return isospinType.at(genus);
}

double ResonanceData::getClebschGordan(int, int , int , int , int , int ) const {
  return 1.;
}

// @TODO Probably take Particles instead of just ids
double ResonanceData::getResonanceSigma(int idA, int idB, double eCM) const {
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  if (eCM < mA + mB) return 0.; // @TODO: This isn't needed if we take two input particles and calculate s

  // Ensure baryon is always idB, if there is a baryon
  if (particleDataPtr->isBaryon(idA) && particleDataPtr->isMeson(idB))
    return getResonanceSigma(idB, idA, eCM);

  auto gens = genify(abs(idA), abs(idB));

  auto totalCharge = particleDataPtr->chargeType(idA) + particleDataPtr->chargeType(idB);
  vector<int> resonanceCandidates;

  if (particleDataPtr->isMeson(idB)) {
    for (auto genus : classToGenera.at("M*")) {
      for (auto species :  genusToParticles.at(genus))
        if (particleDataPtr->chargeType(species) == totalCharge) 
          resonanceCandidates.push_back(species); // break; 
    }
  }
  else {
    for (auto genus : classToGenera.at("N*")) {
      for (auto species :  genusToParticles.at(genus))
        if (particleDataPtr->chargeType(species) == totalCharge) 
          resonanceCandidates.push_back(species); // break; 
    }
    for (auto genus : classToGenera.at("D*")) {
      for (auto species :  genusToParticles.at(genus))
        if (particleDataPtr->chargeType(species) == totalCharge) 
          resonanceCandidates.push_back(species); // break; 
    }
    for (auto genus : classToGenera.at("D")) {
      for (auto species :  genusToParticles.at(genus))
        if (particleDataPtr->chargeType(species) == totalCharge) 
          resonanceCandidates.push_back(species); // break; 
    }

  }

  double sigmaRes = 0;

  cout << " --- Colliding " << particleDataPtr->name(idA) << " with " << particleDataPtr->name(idB) << " --- " << endl;

  for (auto idR : resonanceCandidates) {
    double br = getBR(idR, idA, idB);
    if (br == 0.) continue;

    double gammaRes2 = pow2(particleWidthPtr->mass(speciesToGenus.at(idR), eCM));

    int iA = getIsospin(idA);
    int i3A = 2 * particleDataPtr->chargeType(idA) / 3;
    if (particleDataPtr->isBaryon(idA)) i3A -= 1;
    int iB = getIsospin(idB);
    int i3B = 2 * particleDataPtr->chargeType(idB) / 3;
    if (particleDataPtr->isBaryon(idB)) i3B -= 1;
    int iR = getIsospin(idR);
    int i3R = 2 * particleDataPtr->chargeType(idR) / 3;
    if (particleDataPtr->isBaryon(idR)) i3R -= 1;
    double cg = getClebschGordan(iA, i3A, iB, i3B, iR, i3R);

    int nJRes = particleDataPtr->spinType(idR);
    double m0 = particleDataPtr->m0(idR);

    cout << "Contribution from " << particleDataPtr->name(idR) << endl
          << "Clebsch-Gordan = <" << iA << " " << i3A << " , " << iB << " " << i3B << " | " << iR << " " << i3R << " >^2 " << endl
          << "(2S_R + 1) = " << nJRes << endl
          << "Branching ratio " << br << endl
          << "Gamma = " << sqrt(gammaRes2) << endl
          << endl;

    sigmaRes += pow2(cg) * nJRes * br * gammaRes2/(pow2(eCM - m0) + 0.25 * gammaRes2);
  }

  double s = eCM * eCM;
  double pCMS2 = (s - pow2(mA + mB)) * (s - pow2(mA - mB)) / (4 * s);

  double preFactor = M_PI / (particleDataPtr->spinType(idA) * particleDataPtr->spinType(idB) * pCMS2);
  cout << "Prefactor = " << preFactor << endl;
  sigmaRes *= preFactor;

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
    cout << gen.first << " = { ";
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

void ResonanceData::sanityCheck() {

  for (auto cls : classToGenera) {
    for (auto gen : cls.second) {
      if (genusToParticles.find(gen) == genusToParticles.end())
        cout << " ERROR: particles missing for genus " << gen << " of class " << cls.first << endl;
    }


  }

  for (auto gen : genusToParticles) {
    for (auto spc : gen.second) {
      if (!particleDataPtr->isParticle(spc))
        cout << " ERROR: particle species " << spc << " of genus " << gen.first << " is missing" << endl;
    }
  }

  for (auto clsA : classToGenera) 
  for (auto clsB : classToGenera) {
    pair<ResClass, ResClass> clsPair(clsA.first, clsB.first);
    if (totalSigmaDistribution.find(clsPair) != totalSigmaDistribution.end()) {
      if (scatterChannels.find(clsPair) == scatterChannels.end())
        cout << " ERROR: Total sigma for " << clsPair.first << " " << clsPair.second << " is nonzero, but scattering channel not found" << endl;
      else
      {
        for (auto chn : scatterChannels.at(clsPair))
        for (auto genC : classToGenera.at(chn.first))
        for (auto genD : classToGenera.at(chn.second)) {
          pair<ResGenus, ResGenus> genPair(genC, genD);
          if (partialSigmaDistribution.find(genPair) == partialSigmaDistribution.end())
            cout << " ERROR: Total sigma for " << clsPair.first << " " << clsPair.second << " is nonzero, but relative frequency for output " << genPair.first << " " << genPair.second << " not found" << endl;
        }
      }
    }
  }
}

}
