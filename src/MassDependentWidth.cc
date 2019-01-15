
#include "Pythia8/MassDependentWidth.h"

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

*/
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

bool MassDependentWidth::readXML(istream& stream) {

  string line;

  while (getline(stream, line)) {

    string word1;
    istringstream(line) >> word1;

    if (word1 == "<width") {
      completeTag(stream, line);

      string particleGenus = attributeValue(line, "genus");
      double left = doubleAttributeValue(line, "left");
      double right = doubleAttributeValue(line, "right");

      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);

      this->massDependentWidths.emplace(particleGenus, Interpolator(left, right, data));
    }
    else if (word1 == "<br") {
      completeTag(stream, line);

      string particleGenus = attributeValue(line, "genus");
      string products = attributeValue(line, "products");
      pair<string, string> key(particleGenus, products);
            
      istringstream dataStr(attributeValue(line, "data"));
      vector<double> data;
      double currentData;
      while (dataStr >> currentData)
        data.push_back(currentData);
      
      auto& in = massDependentWidths.at(particleGenus);
      this->branchingRatios.emplace(key, Interpolator(in.left(), in.right(), data));
    }
  }

  return true;
}


double MassDependentWidth::width(string particleGenus, double eCM) const {
  auto entry = massDependentWidths.find(particleGenus);
  if (entry == massDependentWidths.end())
    // @TODO Figure what to do when mass is missing
    return (cout << "Mass-dependent width of particle genus not found: " << particleGenus << endl), 0.;
  else
    return entry->second(eCM);
}

double MassDependentWidth::branchingRatio(string particleGenus, string products, double eCM) const {
  auto entry = branchingRatios.find(pair<string, string>(particleGenus, products));
  if (entry == branchingRatios.end())
    return 0.;
  else
    return entry->second(eCM);
}

const Interpolator& MassDependentWidth::getDistribution(string particleGenus) const {
  return massDependentWidths.at(particleGenus);
}

const Interpolator& MassDependentWidth::getBranchingRatios(string particleGenus, string products) const {
  return branchingRatios.at(pair<string, string>(particleGenus, products));
}

}