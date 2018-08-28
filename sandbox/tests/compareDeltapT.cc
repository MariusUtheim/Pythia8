
#include <iostream>
#include "Pythia8/Pythia.h"
#include "../tests.h"

using namespace Pythia8;

static void readpTs(const Event& evIn, const vector<int> codesIn, 
  vector<int>& countsOut, vector<double>& pTOut) 
{
  for (size_t i = 0; i < evIn.size(); ++i)
  {
    const Particle& p = evIn[i];
    if (!p.isHadron() || !p.isFinal())
      continue;

    for (int i = 0; i < codesIn.size(); ++i)
    {
      if (p.id() != codesIn[i])
        continue;
      
      countsOut[i] += 1;
      pTOut[i] += p.pT();
      break;
    }
  }
}

void test_compareDeltapT()
{
  vector<int> countDef(typeCodes.size());
  vector<double> pTDef(typeCodes.size());

  vector<int> countResc(typeCodes.size());
  vector<double> pTResc(typeCodes.size());


  Pythia pythiaDef;
  pythiaDef.readFile("tests/compareDeltapT.cmnd");
  pythiaDef.readString("Rescattering:rescattering = off");
  pythiaDef.init();
  
  int nEvent = pythiaDef.mode("Main:numberOfEvents");

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythiaDef.next()) continue;
    readpTs(pythiaDef.event, typeCodes, countDef, pTDef);
  }

  Pythia pythiaResc;
  pythiaResc.readFile("tests/compareDeltapT.cmnd");
  pythiaResc.readString("Rescattering:rescattering = on");
  pythiaResc.readString("Rescattering:doSecondRescattering = on");
  pythiaResc.readString("Rescattering:doDecays = off");
  pythiaResc.init();
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythiaResc.next()) continue;
    readpTs(pythiaResc.event, typeCodes, countResc, pTResc);
  }

  for (size_t i = 0; i < typeCodes.size(); ++i) 
  {
    pTDef[i] /= countDef[i]; 
    pTResc[i] /= countResc[i];
  }

  ofstream csvOut("tests/compareDeltapT.csv");

  csvOut << "Type,MC/Data (no rescattering),MC/Data (with rescattering),<pT>(No re),<pT>(Resc),<pT>(data)" << endl;

  for (int i = 0; i < typeCodes.size(); ++i)
    csvOut << typeNames[i] << ","
           << pTDef[i] / datapT[i] << ","
           << pTResc[i] / datapT[i] << ","
           << pTDef[i] << "," << pTResc[i] << "," << datapT[i] << endl;
}