
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
        entry.addProducts(products, br);
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

double ParticleWidths::partialWidth(int id, vector<int> prods, double eCM) const {
  auto iter = entries.find(id);
  return (iter != entries.end()) ? iter->second.getWidth(prods, eCM) : 0.;
}


double ParticleWidths::branchingRatio(int id, vector<int> prods, double eCM) const {
  // @TODO Ordering of products?
  auto iter = entries.find(id);
  return (iter != entries.end()) ? iter->second.getBR(prods, eCM) : 0.;
}

vector<pair<double, vector<int>>> ParticleWidths::getWeightedProducts(int id, double eCM) const {
  auto iter = entries.find(id);
  if (iter == entries.end()) 
    return vector<pair<double, vector<int>>>();
  else {
    const ParticleWidthEntry& entry = iter->second;
    vector<pair<double, vector<int>>> result;
    for (auto prodBRs : entry.branchingRatios) 
      result.push_back(make_pair(prodBRs.second(eCM), prodBRs.first));
    return result;
  }
}

}