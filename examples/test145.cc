#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {

   bool use2to2 = false;
   int nEvent = 100000;
   Pythia pythia;
   Event& event = pythia.event;

   // e+ e- collisions at 1 TeV, stable top, no hadronization
   pythia.readString("Beams:idA =  11");
   pythia.readString("Beams:idB = -11");
   pythia.readString("PDF:lepton = off");
   pythia.readString("Beams:eCM = 1000.");
   if (use2to2) {
     pythia.readString("Top:ffbar2ttbar(s:gmZ) = on");
   } else {
     pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on");
     pythia.readString("32:onMode = off");
     pythia.readString("32:onIfAny = 6");
   }
   pythia.readString("6:mayDecay = off");
   pythia.readString("HadronLevel:Hadronize = off");
   pythia.readString("TimeShower:QEDshowerByQ = off");
   double dead_cone_angle = 170.0/500.0;

   pythia.init();
   Hist dead_cone_ini("Dead Cone Initial Region?", 40, 0.0, 4.0);
   Hist dead_cone_fin("Dead Cone Final Region?", 40, 0.0, 4.0);
   Hist deadIni("Energy flow around initial t", 40, 0.0, 4.0);
   Hist deadFin("Energy flow around final t"  , 40, 0.0, 4.0);
   Hist eInIni("energy inside initial dead cones", 100, 0., 200.);
   Hist eOutIni("energy outside initial dead cones", 100, 0., 200.);
   Hist eInFin("energy inside final dead cones", 100, 0., 200.);
   Hist eOutFin("energy outside final dead cones", 100, 0., 200.);
   Hist nGlu("number of final-state 'gluons'", 100, -0.5, 99.5);
   Hist mGlu("invariant mass of 'gluons'", 40, 0., 400.);

   for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // identify initial and final top quarks; initiate gluon sums
      int iBeg = (use2to2) ? 5 : 6;
      Vec4 ti1 = event[iBeg].p();
      Vec4 ti2 = event[iBeg + 1].p();
      int itf1 = iBeg;
      int itf2 = iBeg;
      for (int i = iBeg + 1; i < event.size(); ++i) {
        if (event[i].id() ==  6) itf1 = i;
        if (event[i].id() == -6) itf2 = i;
      }
      Vec4 tf1 = event[itf1].p();
      Vec4 tf2 = event[itf2].p();
      double eII = 0.;
      double eOI = 0.;
      double eIF = 0.;
      double eOF = 0.;
      Vec4 gSumI, gSumF, gSumA;
      int    nG = 0;

      // find all final particles within angle of 1.0 around top quark
      // and sum into gluon candidate, making sure not to include top itself
      for (int i = 1; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].idAbs() != 6) {
        double thetaI1   = theta(event[i].p(), ti1);
        double thetaI2   = theta(event[i].p(), ti2);
        double thetaIMin = min( thetaI1, thetaI2 ) / dead_cone_angle;
        double thetaF1   = theta(event[i].p(), tf1);
        double thetaF2   = theta(event[i].p(), tf2);
        double thetaFMin = min( thetaF1, thetaF2 ) / dead_cone_angle;
        if (thetaIMin < 2.) deadIni.fill( pow2(thetaIMin), event[i].e() );
        if (thetaFMin < 2.) deadFin.fill( pow2(thetaFMin), event[i].e() );
        if (thetaIMin < 1.) eII += event[i].e();
        else                eOI += event[i].e();
        if (thetaFMin < 1.) eIF += event[i].e();
        else                eOF += event[i].e();
        if (thetaI1 < 1.) gSumI += event[i].p();
        if (thetaF1 < 1.) gSumF += event[i].p();
        ++nG;
        gSumA += event[i].p();
      }

      eInIni.fill( eII );
      eOutIni.fill( eOI );
      eInFin.fill( eIF );
      eOutFin.fill( eOF );
      nGlu.fill( nG );
      if (gSumA.mCalc() > 0.1) mGlu.fill( gSumA.mCalc() );

      // define angle to top direction, define dead cone angle
      // and fill angle squared
      double angleI = theta(ti1,gSumI) / dead_cone_angle;
      if (gSumI.e() != 0) dead_cone_ini.fill( pow2(angleI) );
      double angleF = theta(tf1,gSumF) / dead_cone_angle;
      if (gSumF.e() != 0) dead_cone_fin.fill( pow2(angleF) );
   }

   pythia.stat();
   mGlu /= nEvent;
   cout << dead_cone_ini << dead_cone_fin << deadIni << deadFin
        << eInIni << eOutIni << eInFin << eOutFin << nGlu << mGlu;
  return 0;
}
