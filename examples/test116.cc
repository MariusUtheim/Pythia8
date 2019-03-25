#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  Pythia pythia;
  Event& event = pythia.event;

  //e+ e- inital state
  pythia.readString("Beams:eCM = 1000.");
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");

  //Hidden Valley parameter
  pythia.readString("HiddenValley:Ngauge = 1");
  pythia.readString("HiddenValley:ffbar2Zv = on");
  pythia.readString("HiddenValley:FSR = on");
  //pythia.readString("HiddenValley:alphaFSR = 0.2");
  pythia.readString("HiddenValley:alphaFSR = 0.02");

  //switch off all other kind of radiation
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:Remnants = off");
  pythia.readString("HadronLevel:all = off");
  pythia.readString("Check:event = off");

  // DM mass
  pythia.readString("4900101:m0 = 4.0");

  // Z' mass and decay channel to DM
  pythia.readString("4900023:m0 = 1000.");
  pythia.readString("4900023:onMode = off");
  pythia.readString("4900023:OnIfMatch = 4900101 -4900101");

  // A' mass
  pythia.readString("4900022:m0 = 1.5");
  pythia.readString("4900022:onMode = off");

  pythia.init();

  //histogram
  Hist hist_N("N_A", 100, -0.5, 99.5);
  Hist hist_E("E_A'", 100, 0, 500.);

  //event generation
  double Ngenerate = 400000;
  for (int iEvent = 1; iEvent <= Ngenerate; ++iEvent) {
    if( !pythia.next() ) continue;

    int N_Ap=0; //number A's
    double E_Ap=0; //energy of A'
    double E_DM=0; //energy of highest energetic DM particle

    //loop over particles
    for ( int i = 0; i < event.size(); ++i ) {

      // Count number of final state A's and save one energy
      if( event[i].isFinal() && event[i].id() == 4900022) {
        ++N_Ap;
        E_Ap = event[i].e();
      }

      // Save energy of the hightest energetic dark matter particle.
      //?? if( event[i].isFinal() && event[i].id() == 4900101)
      if( event[i].isFinal() && event[i].idAbs() == 4900101)
        E_DM = (event[i].e() > E_DM) ? event[i].e() : E_DM;
    }

    // Select only events with exactly 1 A' .
    // Including the second condition will reproduce our 'Pythia-A-500' plot.
    // Only first condition reproduces the 'Pythia-A' plot with its high tail.
    if( N_Ap == 1 /* && E_DM > 495. */ ) hist_E.fill(E_Ap);
    hist_N.fill( N_Ap );

  // End loop over events.
  }

  // Statistics and histogram.
  pythia.stat();
  hist_E.takeLog();
  cout << hist_N << hist_E;

  // Done.
  return 0;
}
