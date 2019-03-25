// test243.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Development of the LowEnergyHadHad class for low-energy hadron-hadron
// collisions, as needed e.g. in rescattering descriptions.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// Header file for the LowEnergyHadHad class.
// LowEnergyHadHad: described the low-energy collision between two hadrons.

class LowEnergyHadHad {

public:

  // Constructor.
  LowEnergyHadHad() : infoPtr(), particleDataPtr(), rndmPtr() {}

  // Initialize the class.
  bool init(Info* infoPtrIn, Settings& settings,
    ParticleData* particleDataPtrIn, Rndm* rndmPtrIn);

  // Produce outgoing primary hadrons frrom collision of incoming pair.
  bool collide( int type, int i1, int i2, Event& event);

private:

  // Constants: could only be changed in the code itself.
  static const int MAXLOOP;
  static const double MASSREDUCERATE, MDIFFMIN;

  // Properties of the current collision.
  int    sizeOld;
  double fracetass, fracetaPss, xPowMes, xPowBar, xDiqEnhance, sigmaQ,
         m1, m2, eCM;

  // Pointer to various information on the generation.
  Info*         infoPtr;

  // Pointer to the particle data table.
  ParticleData* particleDataPtr;

  // Pointer to the random number generator.
  Rndm*         rndmPtr;

  // Handle inelastic nondiffractive collision.
  bool nondiff( int i1, int i2, Event& event);

  // Handle elastic and diffractive collisions.
  bool eldiff( int type, int i1, int i2, Event& event);

  // Split a hadron inte a colour and an anticolour part.
  pair< int, int> splitHad( int id);  

  // Choose relative momentum of colour and anticolour constituents in hadron.
  double splitZ( int id1, int id2, double mRat1, double mRat2);

};

//==========================================================================

// The LowEnergyHadHad class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to split beam particles before reconnection.
const int LowEnergyHadHad::MAXLOOP = 100;

// Gradually reduce assumed quark masses from their constituent values.
const double LowEnergyHadHad::MASSREDUCERATE = 0.02;

// Let diffractive mass spectrum begin this much above unexcited mass.
const double LowEnergyHadHad::MDIFFMIN = 0.2;

//--------------------------------------------------------------------------

// Initialize the LowEnergyHadHad class as required.

bool LowEnergyHadHad::init(Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;

  // Mixing for eta and eta'.
  double theta    = settings.parm("StringFlav:thetaPS");
  double alpha    = (theta + 54.7) * M_PI / 180.; 
  fracetass       = pow2(sin(alpha));
  fracetaPss      = 1. - fracetass; 

  // Transverse momentum spread.
  sigmaQ          = settings.parm("StringPT:sigma") / sqrt(2.);

  // Longitudinal momentum sharing of valence quarks in hadrons.
  xPowMes         = settings.parm("BeamRemnants:valencePowerMeson");
  xPowBar         = 0.5 * ( settings.parm("BeamRemnants:valencePowerUinP")
                          + settings.parm("BeamRemnants:valencePowerDinP") );
  xDiqEnhance     = settings.parm("BeamRemnants:valenceDiqEnhance");

  // Done.
  return true;
 
}

//--------------------------------------------------------------------------

// Produce outgoing primary hadrons from collision of incoming pair.
// type = 0: mix; = 1: nondiff; = 2 : el; = 3: SD (XB); = 4: SD (AX); = 5: DD.

bool LowEnergyHadHad::collide( int typeIn, int i1, int i2, Event& event) {

  // Check that incoming hadrons. Store current event size.
  if (!event[i1].isHadron() || !event[i2].isHadron()) return false;
  if (typeIn < 0 || typeIn > 5) return false;
  sizeOld = event.size();

  // Pick event type for typeIn = 0.
  int type = typeIn;
  if (typeIn == 0) {
    // To be done.
  } 

  //  Hadron masses and collision invariant mass.
  m1      = event[i1].m();
  m2      = event[i2].m();
  eCM     = (event[i1].p() + event[i2].p()).mCalc();

  // Do inelastic nondiffractive collision.
  if (type == 1 && !nondiff( i1, i2, event)) return false;

  // Do elastic or diffractive collision.
  if (type > 1 && !eldiff( type, i1, i2, event)) return false;
  
  // Find and do boost from collision rest frame to event frame.
  RotBstMatrix MfromCM = fromCMframe( event[i1].p(), event[i2].p());
  for (int i = sizeOld; i < event.size(); ++i) event[i].rotbst( MfromCM); 

  // Done.
  return true;
 
}

//--------------------------------------------------------------------------

// Do an inelastic nondiffractive scattering.

bool LowEnergyHadHad::nondiff( int i1, int i2, Event& event) {

  // Temporary variables.
  pair< int, int> pair1, pair2;
  int idc1, idac1, idc2, idac2;
  pair<double, double> gauss2;
  double redGrad, redNow, px1, py1, pTs1, px2, py2, pTs2, mc1, mac1, mc2, mac2, 
    mT2c1, mT2ac1, mT2c2, mT2ac2, z1, z2, mT1, mT2; 

  // Check that not stuck in infinite loop.
  int loop = 0;
  do {
    if (++loop == MAXLOOP) {
      infoPtr->errorMsg("Error in LowEnergyHadHad::nondiff: " 
        " failed to construct valid kinematics");  
      return false;
    }
    redGrad = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9)); 

    // Split up hadrons into q + qbar for meson and q + qq for baryon.
    pair1  = splitHad( event[i1].id() );
    idc1   = pair1.first;
    idac1  = pair1.second;
    pair2  = splitHad( event[i2].id() );
    idc2   = pair2.first;
    idac2  = pair2.second;

    // Assign pT kick between constituents inside each hadron.
    gauss2 = rndmPtr->gauss2();
    px1    = sigmaQ * gauss2.first;
    py1    = sigmaQ * gauss2.second; 
    pTs1   = px1 * px1 + py1 * py1;
    gauss2 = rndmPtr->gauss2();
    px2    = sigmaQ * gauss2.first;
    py2    = sigmaQ * gauss2.second; 
    pTs2   = px2 * px2 + py2 * py2;

    // Assign parton masses and squared transverse masses.
    mc1    = particleDataPtr->m0( idc1);
    mac1   = particleDataPtr->m0( idac1);
    redNow = redGrad * min( mA / (mc1 + mac1), 1.);
    mc1   *= redNow;
    mac1  *= redNow;
    mc2    = particleDataPtr->m0( idc2);
    mac2   = particleDataPtr->m0( idac2);
    redNow = redGrad * min( mB / (mc2 + mac2), 1.);
    mc2   *= redNow;
    mac2  *= redNow;
    mT2c1  = pow2(mc1)  + pTs1;
    mT2ac1 = pow2(mac1) + pTs1;
    mT2c2  = pow2(mc2)  + pTs2;
    mT2ac2 = pow2(mac2) + pTs2;

    // Assign relative sharing of longitudinal momentum.
    z1     = splitZ( idc1, idac1, sqrt(mT2c1) / eCM, sqrt(mT2ac1) / eCM); 
    z2     = splitZ( idc2, idac2, sqrt(mT2c2) / eCM, sqrt(mT2ac2) / eCM); 
    mT1    = sqrt( mT2c1 / z1 + mT2ac1 / (1. - z1));
    mT2    = sqrt( mT2c2 / z2 + mT2ac2 / (1. - z2));

  // Ensure that hadron beam remnants are not too massive.
  } while (mT1 + mT2 > eCM);

  // Set up kinematics for outgoing beam remnants. 
  double e1    = 0.5 * (eCM * eCM + mT1 * mT1 - mT2 * mT2) / eCM;
  double pz1   = sqrtpos(e1 * e1 - mT1 * mT1);
  double epz1  = z1 * (e1 + pz1);
  double pzc1  = 0.5 * (epz1 - mT2c1 / epz1 );
  double ec1   = 0.5 * (epz1 + mT2c1 / epz1 ); 
  Vec4 pc1(   px1,  py1,       pzc1,      ec1 );
  Vec4 pac1( -px1, -py1, pz1 - pzc1, e1 - ec1 );    
  double epz2  = z2 * (eCM - e1 + pz1); 
  double pzc2  = -0.5 * (epz2 - mT2c2 / epz2 );
  double ec2   =  0.5 * (epz2 + mT2c2 / epz2 ); 
  Vec4 pc2(   px2,  py2,        pzc2,            ec2 );
  Vec4 pac2( -px2, -py2, -pz1 - pzc2, eCM - e1 - ec2 );    

  // Store new reconnected string systems and mark hadrons decayed.
  int newc = event.nextColTag();  
  event.append(  idc1, 63, i1, 0, 0, 0, newc,     0,  pc1,  mc1);
  event.append( idac2, 63, i2, 0, 0, 0, 0,     newc, pac2, mac2);
  event.append(  idc2, 63, i2, 0, 0, 0, newc + 1, 0,  pc2,  mc2);
  event.append( idac1, 63, i1, 0, 0, 0, 0, newc + 1, pac1, mac1);
  event[i1].statusNeg();
  event[i2].statusNeg();

  // Done.
  return true;
 
}

//--------------------------------------------------------------------------

// Do an elastic or diffractive scattering.
// type = 2: elastic; = 3: SD (XB); = 4: SD (AX); = 5: DD.

bool LowEnergyHadHad::eldiff( int type, int i1, int i2, Event& event) {

  // Classify process type. Find excited mass ranges.
  bool exciteA = (type == 3 || type == 5);
  bool exciteB = (type == 4 || type == 5);
  double mA    = m1;
  double mB    = m2;
  double mAmin = (exciteA) ? m1 + MDIFFMIN : m1;
  double mBmin = (exciteB) ? m2 + MDIFFMIN : m2;
  double mAmax = eCM - mBmin;
  double mBmax = eCM - mAmin;
  if (mAmin + mBmin > eCM) {
    infoPtr->errorMsg("Error in LowEnergyHadHad::eldiff: " 
      " too low invariant mass for diffraction");  
    return false;
  }

  // Temporary variables.
  pair< int, int> pair1, pair2;
  int idc1, idac1, idc2, idac2;
  pair<double, double> gauss2;
  double  redGrad, redNow, px1, py1, pTs1, px2, py2, pTs2, mc1, mac1, 
    mc2, mac2, mT2c1, mT2ac1, mT2c2, mT2ac2, z1, z2, mT1, mT2; 

  // Check that not stuck in infinite loop.
  int loop = 0;
  do {
    if (++loop == MAXLOOP) {
      infoPtr->errorMsg("Error in LowEnergyHadHad::eldiff: " 
        " failed to construct valid kinematics");  
      return false;
    }
    redGrad = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9)); 

    // Split up hadron A into q + qbar for meson and q + qq for baryon.
    // Assign pT kicks and parton (transverse) masses.
    if (exciteA) {
      pair1  = splitHad( event[i1].id() );
      idc1   = pair1.first;
      idac1  = pair1.second;
      gauss2 = rndmPtr->gauss2();
      px1    = sigmaQ * gauss2.first;
      py1    = sigmaQ * gauss2.second; 
      pTs1   = px1 * px1 + py1 * py1;
      mc1    = particleDataPtr->m0( idc1);
      mac1   = particleDataPtr->m0( idac1);
      redNow = redGrad * min( mA / (mc1 + mac1), 1.);
      mc1   *= redNow;
      mac1  *= redNow;
      mc1   *= redfac * massRed;
      mac1  *= redfac * massRed;
      mT2c1  = pow2(mc1)  + pTs1;
      mT2ac1 = pow2(mac1) + pTs1;
    
      // Assign excited mass and check whether kinematics works. 
      mA     = mAmin * pow( mAmax / mAmin, rndmPtr->flat() );
      if (mA < sqrt(mT2c1) + sqrt(mT2ac1)) continue;
    }

    // Split up hadron B into q + qbar for meson and q + qq for baryon.
    // Assign pT kicks and parton (transverse) masses.
    if (exciteB) {
      pair2  = splitHad( event[i2].id() );
      idc2   = pair2.first;
      idac2  = pair2.second;
      gauss2 = rndmPtr->gauss2();
      px2    = sigmaQ * gauss2.first;
      py2    = sigmaQ * gauss2.second; 
      pTs2   = px2 * px2 + py2 * py2;
      mc2    = massred * particleDataPtr->m0( idc2);
      mac2   = massred * particleDataPtr->m0( idac2);  
      redNow = redGrad * min( mB / (mc2 + mac2), 1.);
      mc2   *= redNow;
      mac2  *= redNow;
      mT2c2  = pow2(mc2)  + pTs2;
      mT2ac2 = pow2(mac2) + pTs2;
    
      // Assign excited mass and check whether kinematics works. 
      mB     = mBmin * pow( mBmax / mBmin, rndmPtr->flat() );
      if (mB < sqrt(mT2c2) + sqrt(mT2ac2)) continue;
    } 

  // End loop over tries.
  } 

  // Done.
  return true;
    
}

//-------------------------------------------------------------------------

// Split up a hadron into a colour and an anticolour part, of q or qq kinds.

pair< int, int> LowEnergyHadHad::splitHad( int id) {
   
  // Hadron flavour content.
  int idAbs = abs(id);
  int id1   = (idAbs/1000)%10;
  int id2   = (idAbs/100)%10;
  int id3   = (idAbs/10)%10;
  int id4, id5;

  // Nondiagonal mesons.
  if (id1 == 0 && id2 != id3) {
    if (id != 130 && id != 310) {
      if (id2%2 == 1) swap( id2, id3);
      if (id > 0) return make_pair( id2, -id3);
      else        return make_pair( id3, -id2);

    // K0S and K0L are mixes d sbar and dbar s.   
    } else {
      if (rndmPtr->flat() < 0.5) return make_pair( 3, -1);
      else                       return make_pair( 1, -3);
    }
      
  // Diagonal mesons: assume complete mixing ddbar and uubar.
  } else if (id1 == 0) {
   if (id2 < 3 || id == 331) {
     id4 = (rndmPtr->flat() < 0.5) ? 1 : 2; 
     // eta and eta' can also be s sbar.
     if (id == 221 && rndmPtr->flat() < fracetass) id4 = 3;
     if (id == 331 && rndmPtr->flat() < fracetaPss) id4 = 3;
     return make_pair( id4, -id4);
   }

  // Octet baryons.
  } else if (idAbs%10 == 2) {
    // Three identical quarks: emergency in case of higher spin 1/2 multiplet.
    if (id1 == id2 && id2 == id3) {id4 = id1; id5 = 1100 * id1 + 3;}
    // Two identical quarks, like normal p or n.
    else if (id1 == id2 || id2 == id3) {
      double rr6 = 6. * rndmPtr->flat();
      if    (id1 == id2 && rr6 < 2.) { id4 = id3; id5 = 1100 * id1 + 3;}
      else if             (rr6 < 2.) { id4 = id1; id5 = 1100 * id3 + 3;}
      else if (rr6 < 3.) { id4 = id2; id5 = 1000 * id1 + 100 * id3 + 3;}
      else               { id4 = id2; id5 = 1000 * id1 + 100 * id3 + 1;}
    // Three nonidentical quarks, Sigma- or Lambda-like.
    } else { 
      int isp = (id2 > id3) ? 3 : 1;    
      if (id3 > id2) swap( id2, id3);   
      double rr12 = 12. * rndmPtr->flat();
      if      (rr12 < 4.) { id4 = id1; id5 = 1000 * id2 + 100 * id3 + isp;}   
      else if (rr12 < 5.) { id4 = id2; id5 = 1000 * id1 + 100 * id3 + isp;}
      else if (rr12 < 6.) { id4 = id3; id5 = 1000 * id1 + 100 * id2 + isp;}
      else if (rr12 < 9.) { id4 = id2; id5 = 1000 * id1 + 100 * id3 + 4 - isp;}
      else                { id4 = id3; id5 = 1000 * id1 + 100 * id2 + 4 - isp;}
    }
    if (id > 0) return make_pair(  id4,  id5);
    else        return make_pair( -id5, -id4);
    
  // Decuplet baryons.
  } else {
    double rr3 = 3. * rndmPtr->flat();
    if (rr3 < 1.)      { id4 = id1; id5 = 1000 * id2 + 100 * id3 + 3;}
    else if (rr3 < 2.) { id4 = id2; id5 = 1000 * id1 + 100 * id3 + 3;} 
    else               { id4 = id3; id5 = 1000 * id1 + 100 * id2 + 3;} 
    if (id > 0) return make_pair(  id4,  id5);
    else        return make_pair( -id5, -id4); 
  }

  // Done. (Fake call to avoid unwarranted compiler warning.)
  return make_pair( 0, 0);

}

//-------------------------------------------------------------------------

// Find relative momentum of colour and anticolour constituents in hadron.

double LowEnergyHadHad::splitZ( int id1, int id2, double mRat1, double mRat2) {

  // Initial values.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  if (id2Abs > 10) swap( mRat1, mRat2);
  double x1, x2, x1a, x1b;

  // Handle mesons.
  if (id1Abs < 10 && id1Abs < 10) {
    do x1 = pow2( mRat1 + (1. - mRat1) * rndmPtr->flat() );
    while ( pow(1. - x1, xPowMes) < rndmPtr->flat() );
    do x2 = pow2( mRat2 + (1. - mRat2) * rndmPtr->flat() );
    while ( pow(1. - x2, xPowMes) < rndmPtr->flat() );

  // Handle baryons.
  } else {
    double mRat1ab = 0.5 * mRat1 / xDiqEnhance;
    do x1a = pow2( mRat1ab + (1. - mRat1ab) * rndmPtr->flat() );
    while ( pow(1. - x1a, xPowBar) < rndmPtr->flat() );
    do x1b = pow2( mRat1ab + (1. - mRat1ab) * rndmPtr->flat() );
    while ( pow(1. - x1b, xPowBar) < rndmPtr->flat() );
    x1 = xDiqEnhance * ( x1a + x1b);
    do x2 = pow2( mRat2 + (1. - mRat2) * rndmPtr->flat() );
    while ( pow(1. - x2, xPowBar) < rndmPtr->flat() );
    if (id2Abs > 10) swap( x1, x2);
  }

  // Return z value.      
  return x1 / (x1 + x2); 
  
}

//==========================================================================

int main() {

  // Number of events to generate. Max number of errors.
  int nEvent = 2;
  int nAbort = 5;

  // Set up inelastic nondiffractive LHC collisions.
  Pythia pythia;
  Event& event = pythia.event;
  pythia.readString("SoftQCD:nonDiffractive = on");
  pythia.init();

  // Initialize LowEnergyHadHad class.
  LowEnergyHadHad LEHH;
  LEHH.init( &pythia.info, pythia.settings, &pythia.particleData,
    &pythia.rndm); 

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

   // Find a hadron pair with invariant mass above 10 GeV.
   int ncollide = 0;
   for (int i1 = 0; i1 < event.size() - 1; ++i1) 
   if (event[i1].isFinal() && event[i1].isHadron()) {
     for (int i2 = i1 + 1; i2 < event.size(); ++i2) 
     if (event[i2].isFinal() && event[i2].isHadron() 
     && (event[i1].p() + event[i2].p()).mCalc() > 10.) {
       cout << " collision between " << i1 << " and " << i2 << endl;
       LEHH.collide( i1, i2, event);
       ++ncollide;
       break;
     }
   } 
   cout << " event contained " << ncollide << " collisions " << endl;
   event.list();

  // End event loop. Final statistics.
  }
  pythia.stat();

  // Done.
  return 0;
}
