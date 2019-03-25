#include "Pythia8/Pythia.h"
#include <cassert>
using namespace Pythia8;

//==========================================================================

// Write own derived UserHooks class.

class MyUserHooks : public UserHooks {

public:

  // Constructor and destructor.
  MyUserHooks(bool vetoG2GGIn, bool vetoG2QQIn) {vetoG2GG = vetoG2GGIn;
    vetoG2QQ = vetoG2QQIn; }
  ~MyUserHooks() {}

  // Allow FSR emissions to be vetoed.
  virtual bool canVetoFSREmission() {return true;}

  // Veto emissions if vetoing is on and branching parton is a gluon,
  // either to gluons or to quarks.
  virtual bool doVetoFSREmission( int sizeOld, const Event& event, int,
    bool = false) {
    if (event[event[sizeOld].mother1()].id() == 21) {
      if (vetoG2GG && event[sizeOld].id() == 21) return true;
      if (vetoG2QQ && event[sizeOld].idAbs() < 6) return true;
    }
    return false;
  }

private:
  bool vetoG2GG, vetoG2QQ;

};

//==========================================================================

// Method to switch distance measure between e+e- and pp.

double sep( Vec4 p1, Vec4 p2, bool useEE) {
  if (useEE) return theta( p1, p2);
  else       return RRapPhi( p1, p2);
}

//==========================================================================

int main() {

   // Main setup choices.
   // useLC : e+e- linear collider, else LHC.
   bool useLC      = false;
   // use2to2: true if LC ttbar set up as 2 -> 2, else as 2 -> 1 -> 2.
   bool use2to2    = false;
   // vetoGluons: forbid gluons from branching to gluons or quarks.
   bool vetoG2GG   = false;
   bool vetoG2QQ   = false;
   // deadCone: use eikonal to suppress radiation in dead cone.
   bool deadCone   = true;
   // meForTT: find ME for ttbar pair no matter what.
   bool meForTT    = false;

   // Number of events.
   int  nEvent     = 1000000;

   // Generator. User hook to handle vetoing. Dead-cone setting.
   Pythia pythia;
   Event& event = pythia.event;
   MyUserHooks* myUserHooks = new MyUserHooks( vetoG2GG, vetoG2QQ);
   pythia.setUserHooksPtr( myUserHooks);
   pythia.settings.flag("TimeShower:recoilDeadCone", deadCone);
   pythia.settings.flag("TimeShower:MEextended", meForTT);

   // e+ e- collisions at 1 TeV.
   if (useLC) {
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

   // pp collisions at 13 TeV LHC.
   } else {
     pythia.readString("Beams:idA = 2212");
     pythia.readString("Beams:idB = 2212");
     pythia.readString("Beams:eCM = 13000.");
     pythia.readString("Top:gg2ttbar = on");
     pythia.readString("Top:qqbar2ttbar = on");
     pythia.readString("PhaseSpace:pTHatMin = 500.");
   }

   // Stable top, no MPI, no hadronization, no QED shower.
   pythia.readString("6:mayDecay = off");
   pythia.readString("PartonLevel:MPI = off");
   pythia.readString("HadronLevel:Hadronize = off");
   pythia.readString("TimeShower:QEDshowerByQ = off");
   pythia.readString("Next:numberCount = 100000");

   // Special options for tests.
   //pythia.readString(" TimeShower:MEcorrections = off");
   //pythia.readString("TimeShower:alphaSorder = 0");
   //pythia.readString("TimeShower:alphaSvalue = 0.2");
   //pythia.readString("TimeShower:nGluonToQuark = 0");
   pythia.readString("PartonLevel:ISR = off");

   // Initialize.
   pythia.init();
   double dead_cone_angle = 170.0/500.0;

   // Histograms.
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
   Hist angle_ratio_1("Angle Ratio 1", 50, 0.0, 4.0);
   Hist angle_ratio_2("Angle Ratio 2", 50, 0.0, 4.0);

   // Loop over event generation.
   for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Identify initial and final top quarks; initiate gluon sums
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

      // Find all final particles within angle of 1.0 around top quark
      // and sum into gluon candidate, making sure not to include top itself
      for (int i = 1; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].idAbs() != 6) {
        double thetaI1   = sep(event[i].p(), ti1, useLC);
        double thetaI2   = sep(event[i].p(), ti2, useLC);
        double thetaIMin = min( thetaI1, thetaI2 ) / dead_cone_angle;
        double thetaF1   = sep(event[i].p(), tf1, useLC);
        double thetaF2   = sep(event[i].p(), tf2, useLC);
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
      double angleI = sep(ti1,gSumI,useLC) / dead_cone_angle;
      if (gSumI.e() != 0) dead_cone_ini.fill( pow2(angleI) );
      double angleF = sep(tf1,gSumF,useLC) / dead_cone_angle;
      if (gSumF.e() != 0) dead_cone_fin.fill( pow2(angleF) );

      // Define particles of interest
      Particle initial_top = use2to2 ? event[5] : event[6];
      Particle radiated_gluon;
      Particle split_gluon_1;
      Particle split_gluon_2;

      // first find radiated gluon from top quark
      vector<int> daughters = initial_top.daughterList();
      if (daughters.size() == 0) continue;
      while (daughters.size() == 1) {
         initial_top = event[daughters[0]];
         daughters = initial_top.daughterList();
      }
      if (daughters.size() == 0) continue;
      assert(daughters.size() == 2);
      initial_top = event[daughters[0]];
      radiated_gluon = event[daughters[1]];
      //assert(initial_top.id() == 6);
      assert(initial_top.idAbs() == 6);
      assert(radiated_gluon.id() == 21 || radiated_gluon.id() == 22);

      // next find gluon splitting
      daughters = radiated_gluon.daughterList();
      if (daughters.size() == 0) continue;
      while (daughters.size() == 1) {
         radiated_gluon = event[daughters[0]];
         daughters = radiated_gluon.daughterList();
      }
      if (daughters.size() == 0) continue;
      assert(daughters.size() == 2);
      split_gluon_1 = event[daughters[0]];
      split_gluon_2 = event[daughters[1]];

      // define angle to initial gluon and angle to split gluons
      double angt0 = (sep(initial_top.p(),radiated_gluon.p(),useLC))
                   /(dead_cone_angle);
      double angt1 = (sep(initial_top.p(),split_gluon_1.p(),useLC))
                   /(dead_cone_angle);
      double angt2 = (sep(initial_top.p(),split_gluon_2.p(),useLC))
                   /(dead_cone_angle);
      angle_ratio_1.fill(angt1*angt1/angt0/angt0);
      angle_ratio_2.fill(angt2*angt2/angt0/angt0);
   }

   // Statistics and histograms.
   pythia.stat();
   mGlu /= nEvent;
   cout << dead_cone_ini << dead_cone_fin << deadIni << deadFin
        << eInIni << eOutIni << eInFin << eOutFin << nGlu << mGlu
        << angle_ratio_1 << angle_ratio_2;

  return 0;
}
