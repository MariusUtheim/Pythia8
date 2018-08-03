// ParticleData.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// DecayChannel, ParticleDataEntry and ParticleData classes.

#include "Pythia8/CrossSectionData.h"


namespace Pythia8 {

//==========================================================================

// CrossSectionDataEntry class.
// @TODO comment

//--------------------------------------------------------------------------
  
InteractionChannel& CrossSectionDataEntry::pickChannel() {

  // Find channel in table.
  double brSum = 0;
  for (size_t i = 0; i < channels.size(); ++i)
    if (channels[i].onMode())
      brSum += channels[i].bRatio();

  double brThreshold = brSum * crossSectionDataPtr->rndmPtr->flat();

  for (size_t i = 0; i < channels.size(); ++i) 
    if (channels[i].onMode()) {
      brThreshold -= channels[i].bRatio();
      if (brThreshold <= 0)
        return channels[i];
    }

  // @TODO: This can never be reached (as long as branching ratios are all nonnegative)
  throw "pickChannel() failed";
}

//==========================================================================

// CrossSectionData class.
// @TODO comment

//--------------------------------------------------------------------------

// @TODO comment

void CrossSectionData::initCommon() {
	// @TODO fill or delete
}

//--------------------------------------------------------------------------

// Read in database from specific XML file (which may refer to others).

bool CrossSectionData::readXML(string inFile, bool reset) {
  ifstream is(inFile.c_str());
  return readXML(is, reset);
}

//--------------------------------------------------------------------------

// Read in database from specific XML stream (which may refer to others).

bool CrossSectionData::readXML(istream &inStr, bool reset) {

  // Normally reset whole database before beginning.
  if (reset) {
    pdt.clear();
		// @TODO: further reinitialisation
    isInit = false;
  }

  // Check that instream is OK.
  if (!inStr.good()) {
    infoPtr->errorMsg("Error in CrossSectionData::readXML:"
      " did not find data");
    return false;
  }

  // Read in one line at a time.
  string currentLine;
  vector<string> lines;
  while (getline(inStr, currentLine)) {
  // @TODO: Check that this part is necessary
    // Get first word of a line.
    //istringstream getfirst(currentLine);

    lines.push_back(currentLine);
  }


  // Number of lines saved.
  int nLines = lines.size();

  CrossSectionDataEntry* currentEntry = 0;

	for (int i = 0; i < nLines; ++i) {

    // Retrieve line.
    string line = lines[i];

    // Get first word of a line.
    istringstream getfirst(line);
    string word1;
    getfirst >> word1;

    // Check for occurence of a cross-section data. Add any continuation lines.
    if (word1 == "<crossSection") {
      while (line.find(">") == string::npos) {
        if (++i >= nLines) break;
        string addLine = lines[i];
        line += addLine;
      }

      // Read in particle properties.
			int idATmp = intAttributeValue(line, "idA");
			int idBTmp = intAttributeValue(line, "idB");
			double sigmaTmp = doubleAttributeValue(line, "sigma");

      // Erase if particle already exists.
      if (isCrossSection(idATmp, idBTmp)) 
				pdt.erase(pairify(idATmp, idBTmp));

      // Create new cross section entry
      currentEntry = &addCrossSection(idATmp, idBTmp, sigmaTmp);
      
      // Check for occurence of a decay channel. Add any continuation lines.
    } else if (word1 == "<interactionChannel") {
      while (line.find(">") == string::npos) {
        if (++i >= nLines) break;
        string addLine = lines[i];
        line += addLine;
      } 

      // Read in channel properties - products so far only as a string.
      int resonanceTmp   = intAttributeValue(line, "resonance");
      int onModeTmp      = intAttributeValue(line, "onMode");
      double bRatioTmp   = doubleAttributeValue(line, "bRatio");
      int meModeTmp      = intAttributeValue(line, "meMode");


      if (resonanceTmp == 0) {
        infoPtr->errorMsg("Error in CrossSectionData::readXML:"
                          " no resonance", line);
        return false;
      }

      // Store new channel (if particle already known).
      if (!currentEntry) { // @TODO: This is sensitive to node naming conflicts (e.g. channels from ParticleData)
        infoPtr->errorMsg("Error in CrossSectionData::readXML:"
                          " orphan decay channel", line);
        return false;
      }

      currentEntry->addChannel(resonanceTmp, onModeTmp, bRatioTmp, meModeTmp);
    };
  };

  // All particle data at this stage defines baseline original.
  if (reset) 
    for (map<pair<int, int>, CrossSectionDataEntry>::iterator pdtEntry 
         = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    pdtEntry->second.hasChanged(false);
  }

  // Done.
  isInit = true;
  return true;
}

//--------------------------------------------------------------------------

// Print out complete or changed table of database in numerical order.
void CrossSectionData::list(bool changedOnly, bool changedRes) {
  // @TODO Print the desired format
  // Table header; output for bool as off/on.
  if (!changedOnly) {
    cout << "\n -------  PYTHIA Cross-Section Data Table (complete)  --------"
         << "------------------------------------------------------------"
         << "--------------\n \n";

  } else {
    cout << "\n --------  PYTHIA Cross-Section Data Table (changed only)  ----"
         << "------------------------------------------------------------"
         << "--------------\n \n";
  }
  cout << "      idA  idB           sigma \n"
       << "      | no onMode    bRatio   meMode      products \n";
  
  // Iterate through the particle data table. Option to skip unchanged.
  int nList = 0;
  for (map<pair<int, int>, CrossSectionDataEntry>::iterator pdtIter // @TODO: Iterator name
    = pdt.begin(); pdtIter != pdt.end(); ++pdtIter) {
    CrossSectionDataEntry& entry = pdtIter->second;
    if ( !changedOnly || entry.hasChanged() ||
      ( changedRes && entry.resonancePtr() != 0 ) ) {

      // Print particle properties.
      ++nList;

      cout << "\n" 
        << setw(9) << entry.idA() 
        << setw(5) << entry.idB()
        << setw(16) << scientific << setprecision(5) << entry.sigma()
        << endl;
      
      // Loop through the decay channel table for each particle.
      if (entry.sizeChannels() > 0) { // @TODO: If not, this is an error! And either way, the check is not needed
        for (int i = 0; i < int(entry.sizeChannels()); ++i) {
          InteractionChannel& channel = entry.channel(i);
          cout << "      | "
               << setw(2)  << i
               << setw(7)  << channel.onMode()
               << setw(11) << fixed << channel.bRatio()
               << setw(8)  << channel.meMode() << " "; 

          // @TODO: Fix
          //for (int j = 0; j < channel.nProducts(); ++j)
            //cout << setw(8) << channel.product(j) << " ";

          cout << "\n";
        }
      }
    }

  }

  // End of loop over database contents.
  if (changedOnly && nList == 0) cout << "\n no particle data has been "
       << "changed from its default value \n";
  cout << "\n --------  End PYTHIA Particle Data Table  -----------------"
       << "--------------------------------------------------------------"
       << "----------\n" << endl;

}


//--------------------------------------------------------------------------

// Extract XML value string following XML attribute.

string CrossSectionData::attributeValue(string line, string attribute) {

  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute);
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

}

//--------------------------------------------------------------------------

// Extract XML bool value following XML attribute.

bool CrossSectionData::boolAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return false;
  return boolString(valString);
}

//--------------------------------------------------------------------------

// Extract XML int value following XML attribute.

int CrossSectionData::intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0;
  istringstream valStream(valString);
  int intVal;
  valStream >> intVal;
  return intVal;

}

//--------------------------------------------------------------------------

// Extract XML double value following XML attribute.

double CrossSectionData::doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.;
  istringstream valStream(valString);
  double doubleVal;
  valStream >> doubleVal;
  return doubleVal;

}

//==========================================================================


} // end namespace Pythia8
