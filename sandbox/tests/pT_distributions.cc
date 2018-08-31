
#include "../tests.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

void test_pT_distributions()
{
  vector<int> countsDef(typeCodes.size()), countsResc(typeCodes.size());
  vector<double> pTsDef(typeCodes.size()), pTsResc(typeCodes.size());

  Pythia pythiaDef;
  pythiaDef.readFile("tests/common.cmnd");
  pythiaDef.readString("HadronLevel:Rescatter = off");
  pythiaDef.init();

  int nEvent = pythiaDef.mode("Main:numberOfEvents");

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythiaDef.next()) continue;

    Event& ev = pythiaDef.event;
    for (size_t iParticle = 0; iParticle < ev.size(); ++iParticle)
    {
      Particle&p = ev[iParticle];
      if (!p.isHadron() || !p.isFinal() || abs(p.y()) > 0.5) continue;

      for (int i = 0; i < typeCodes.size(); ++i)
      {
        if (p.id() != typeCodes[i]) continue;

        countsDef[i] += 1;
        pTsDef[i] += p.pT();
        break;
      }
    }
  }

  Pythia pythiaResc;
  pythiaResc.readFile("tests/common.cmnd");
  pythiaResc.readString("HadronLevel:Rescatter = on");
  pythiaResc.init();

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythiaResc.next()) continue;

    Event& ev = pythiaResc.event;
    for (size_t iParticle = 0; iParticle < ev.size(); ++iParticle)
    {
      Particle&p = ev[iParticle];
      if (!p.isHadron() || !p.isFinal() ) continue;

      for (int i = 0; i < typeCodes.size(); ++i)
      {
        if (p.id() != typeCodes[i]) continue;

        countsResc[i] += 1;
        pTsResc[i] += p.pT();
        break;
      }
    }
  }

  ofstream csvOut("tests/pT_distributions.csv");
  csvOut << "Type,<pT>(data),<pT>(Pythia),<pT>(Rescattering)" << endl;

  for (size_t i = 0; i < typeCodes.size(); ++i)
    csvOut << typeNames[i] << "," << datapT[i] << "," 
           << (pTsDef[i] / countsDef[i]) << ","
           << (pTsResc[i] / countsResc[i]) << endl;
}
