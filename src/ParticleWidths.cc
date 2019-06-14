
#include "Pythia8/ParticleWidths.h"

namespace Pythia8 {


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


bool ParticleWidths::readXML(istream& stream) {

  string line;

  while (getline(stream, line)) {

    string word1;
    istringstream(line) >> word1;

    if (word1 == "<width") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");
      double left = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      this->entries.emplace(id, ParticleWidthEntry(Interpolator(left, right, data)));
    }
    else if (word1 == "<br") {
      completeTag(stream, line);

      int id = intAttributeValue(line, "id");

      istringstream productStr(attributeValue(line, "products"));
      vector<int> products;
      int currentProduct;
      while (productStr >> currentProduct)
        products.push_back(currentProduct);

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);
      
      auto iter = entries.find(id);
      if (iter == entries.end()) {
        infoPtr->errorMsg( "Warning in ParticleWidths::readXML: "
          "got br for a particle with undefined width");
        return false;
      }
      else {
        ParticleWidthEntry& entry = iter->second;
        Interpolator br(entry.widths.left(), entry.widths.right(), data);
        // @TODO: Don't stream to a vector; get it directly into a pair
        entry.addProducts(make_pair(products[0], products[1]), br);
      }
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
  auto iter = entries.find(id);
  return (iter != entries.end()) ? iter->second.widths(eCM) : 0.;
}

double ParticleWidths::partialWidth(int id, pair<int, int> prods, double eCM) const {
  auto iter = entries.find(id);
  return (iter != entries.end()) ? iter->second.getWidth(prods, eCM) : 0.;
}


double ParticleWidths::branchingRatio(int id, pair<int, int> prods, double eCM) const {
  // @TODO Ordering of products?
  auto iter = entries.find(id);
  return (iter != entries.end()) ? iter->second.getBR(prods, eCM) : 0.;
}

vector<pair<double, pair<int, int>>> ParticleWidths::getWeightedProducts(int id, double eCM) const {
  auto iter = entries.find(id);
  if (iter == entries.end()) 
    return vector<pair<double, pair<int, int>>>();
  else {
    const ParticleWidthEntry& entry = iter->second;
    vector<pair<double, pair<int, int>>> result;
    for (auto prodBRs : entry.branchingRatios) 
      result.push_back(make_pair(prodBRs.second(eCM), prodBRs.first));
    return result;
  }
}

pair<int, int> ParticleWidths::pickDecayChannel(int idRes, double eCM) {

  auto brs = getWeightedProducts(idRes, eCM);
  if (brs.size() == 0) {
    // @TODO: This would be a bug
    cout << "Got no decay modes" << endl;
    return make_pair(0, 0);
  }

  vector<double> weights(brs.size());
  for (size_t i = 0; i < brs.size(); ++i)
    weights[i] = brs[i].first;

  return brs[rndmPtr->pick(weights)].second;
}


static double pCMS(double eCM, double mA, double mB) {
  double sCM = eCM * eCM;
  return sqrt((sCM - pow2(mA + mB)) * (sCM - pow2(mA - mB))) / (2. * eCM);
}

static double breitWigner(double gamma, double dm) {
  return 1. / (2. * M_PI) * gamma / (dm * dm + 0.25 * gamma * gamma);
}

static constexpr double MAX_LOOPS = 100;

double ParticleWidths::pickMass(int idRes, double eCM, double mB, int lType) {
  auto iter = entries.find(idRes);
  if (iter == entries.end())
    return 0.;

  ParticleWidthEntry& channel(iter->second);

  // @TODO: Maybe an mPeak that is different from m0 will be more efficient
  double mMin = channel.widths.left(), 
         mMax = min(channel.widths.right(), eCM - mB),
         mPeak = channel.m0,
         m0 = channel.m0;
  double gamma = channel.widths(mPeak);

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
                  * breitWigner(channel.widths(mCand), mCand - m0);
  
    if (rndmPtr->flat() * envelope < yDistr)
      return mCand;
  }

  //// Alternative implementation using lambda functions
  //auto sampler = [&]() {
  //  if (rndmPtr->pick(ps) == 0) {
  //    double r = (0.5 - cdfLow) + rndmPtr->flat() * cdfLow;
  //    return mPeak + 0.5 * gamma * tan(M_PI * (r - 0.5));
  //  }
  //  else {
  //    double r = 0.5 + rndmPtr->flat() * cdfHigh;
  //    return mPeak + gamma * tan(M_PI * (r - 0.5));
  //  }
  //};
//
  //auto accepter = [&](double m) {
  //  double overestimate = m < m0 ? scale * breitWigner(gamma, m - mPeak)
  //                      : scale * 2. * breitWigner(2. * gamma, m - mPeak);
  //  double distr = pow(pCMS(eCM, m, mB), lType) 
  //             * breitWigner(channel.widths(m), m - m0);
  //  return distr / overestimate;
  //};
  //if (!rndmPtr->hitAndMiss(distribution, accepter, mCand)) {
  //  infoPtr->errorMsg("Warning in ParticleWidths::pickMass: "
  //    "Could not choose mass within the prescribed number of iterations");
  //  return m0;
  //}
  //
  //return mCand;

  infoPtr->errorMsg("Warning in ParticleWidths::pickMass: "
    "Could not choose mass within the prescribed number of iterations");
  return m0;
}

// @TODO: Implement this properly
//        For the first iteration, we just pick one particle on shell
pair<double, double> ParticleWidths::pickMass2(int id1, int id2, double eCM, int lType) {
  if (id1 > id2)
    swap(id1, id2);

  auto iter = entries.find(id2);
  if (iter == entries.end()) return make_pair(0., 0.);

  double m2 = iter->second.m0;
  double m1 = pickMass(id1, eCM, m2, lType);

  return make_pair(m1, m2);
}

}