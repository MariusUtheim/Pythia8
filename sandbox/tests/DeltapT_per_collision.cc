
#include "../tests.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

void test_DeltapT_per_collision()
{
  vector<int> typeCodes =    { 221,   321,  313,   2212, 333,   3312,   3224,      3324,    };
  vector<string> typeNames = { "pi+", "K+", "K*0", "p",  "phi", "Chi-", "Sigma*+", "Chi*0", };

  vector<int> counts(typeCodes.size());
  vector<double> pTs(typeCodes.size());

  Pythia pythia;
  pythia.readFile("tests/DeltapT_per_collision.cmnd");
  pythia.init();

  int nEvent = pythia.mode("Main:numberOfEvents");

  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythia.next()) continue;

    Event& ev = pythia.event;
    for (size_t iParticle = 0; iParticle < ev.size(); ++iParticle)
    {
      Particle&p = ev[iParticle];
      if (p.isFinal()
        || !(ev[p.daughter1()].statusAbs() == 111 || ev[p.daughter1()].statusAbs() == 112)
        || (ev[p.daughter1()].id() == ev[p.daughter2()].id())
      ) continue;
      
      for (size_t i = 0; i < typeCodes.size(); ++i)
      {
        if (p.id() != typeCodes[i])
          continue;
        int daughterId = ev[p.daughter1()].id() == p.id()
                          ? p.daughter1() 
                          : p.daughter2();

        double pT0 = p.pT(), pT1 = ev[daughterId].pT();
        
        counts[i] += 1;
        pTs[i] += (pT1 - pT0);

        break;
      }
    }
  }

  ofstream csvOut("tests/DeltapT_per_collision.csv");
  csvOut << "Type,DeltapT" << endl;

  for (size_t i = 0; i < typeCodes.size(); ++i)
    csvOut << typeNames[i] << "," << (pTs[i] / counts[i]) << endl;
}