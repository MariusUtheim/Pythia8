// test251.cc is a part of the PYTHIA event generator.
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
  int    sizeOld, id1, id2, idc1, idac1, idc2, idac2;

  double fracEtass, fracEtaPss, xPowMes, xPowBar, xDiqEnhance, sigmaQ,
         m1, m2, eCM, sCM, z1, z2, mT1, mT2, mA, mB, 
         mc1, mac1, px1, py1, pTs1, mTsc1, mTsac1, mTc1, mTac1,
         mc2, mac2, px2, py2, pTs2, mTsc2, mTsac2, mTc2, mTac2;

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

  // Split up hadron A or B into a colour pair, with masses and pT values.
  bool splitA( double redMpT);
  bool splitB( double redMpT); 

  // Split a hadron inte a colour and an anticolour part.
  pair< int, int> splitFlav( int id);  

  // Choose relative momentum of colour and anticolour constituents in hadron.
  double splitZ( int iq1, int iq2, double mRat1, double mRat2);

  // Pick slope b of exp(b * t) for elastic and diffractive events.
  double bSlope();

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
  fracEtass       = pow2(sin(alpha));
  fracEtaPss      = 1. - fracEtass; 

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
  id1     = event[i1].id();
  id2     = event[i2].id();
  m1      = event[i1].m();
  m2      = event[i2].m();
  eCM     = (event[i1].p() + event[i2].p()).mCalc();
  sCM     = eCM * eCM;

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

  // Check that not stuck in infinite loop.
  int loop = 0;
  do {
    if (++loop == MAXLOOP) {
      infoPtr->errorMsg("Error in LowEnergyHadHad::nondiff: " 
        " failed to construct valid kinematics");  
      return false;
    }
    double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9)); 

    // Split up hadrons A  and B into q + qbar for meson and q + qq for baryon.
    splitA( redStep); 
    splitB( redStep); 

    // Assign relative sharing of longitudinal momentum.
    z1     = splitZ( idc1, idac1, mTc1 / eCM, mTac1 / eCM); 
    z2     = splitZ( idc2, idac2, mTc2 / eCM, mTac2 / eCM); 
    mT1    = sqrt( mTsc1 / z1 + mTsac1 / (1. - z1));
    mT2    = sqrt( mTsc2 / z2 + mTsac2 / (1. - z2));

  // Ensure that hadron beam remnants are not too massive.
  } while (mT1 + mT2 > eCM);

  // Set up kinematics for outgoing beam remnants. 
  double e1    = 0.5 * (sCM + mT1 * mT1 - mT2 * mT2) / eCM;
  double pz1   = sqrtpos(e1 * e1 - mT1 * mT1);
  double epz1  = z1 * (e1 + pz1);
  double pzc1  = 0.5 * (epz1 - mTsc1 / epz1 );
  double ec1   = 0.5 * (epz1 + mTsc1 / epz1 ); 
  Vec4 pc1(   px1,  py1,       pzc1,      ec1 );
  Vec4 pac1( -px1, -py1, pz1 - pzc1, e1 - ec1 );    
  double epz2  = z2 * (eCM - e1 + pz1); 
  double pzc2  = -0.5 * (epz2 - mTsc2 / epz2 );
  double ec2   =  0.5 * (epz2 + mTsc2 / epz2 ); 
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
  bool excite1 = (type == 3 || type == 5);
  bool excite2 = (type == 4 || type == 5);
  mA           = m1;
  mB           = m2;
  double mAmin = (excite1) ? m1 + MDIFFMIN : m1;
  double mBmin = (excite2) ? m2 + MDIFFMIN : m2;
  double mAmax = eCM - mBmin;
  double mBmax = eCM - mAmin;
  if (mAmin + mBmin > eCM) {
    infoPtr->errorMsg("Error in LowEnergyHadHad::eldiff: " 
      " too low invariant mass for diffraction");  
    return false;
  }

  // Check that not stuck in infinite loop.
  int loop = 0;
  do {
    if (++loop == MAXLOOP) {
      infoPtr->errorMsg("Error in LowEnergyHadHad::eldiff: " 
        " failed to construct valid kinematics");  
      return false;
    }
    double redStep = (loop < 10) ? 1. : exp( -MASSREDUCERATE * (loop - 9)); 

    // Split up hadron 1 to 3 (side A) and assign excited mass.
    if (excite1) {
      splitA( redStep); 
      mA     = mAmin * pow( mAmax / mAmin, rndmPtr->flat() );
      if (mA < mTc1 + mTac1) continue;
    }

    // Split up hadron 2 to 4 (side B) and assign excited mass.
    if (excite2) {
      splitB( redStep);
      mB     = mBmin * pow( mBmax / mBmin, rndmPtr->flat() );
      if (mB < mTc2 + mTac2) continue;
    } 

  // Ensure that pair of hadron masses not too large. Squared masses.
  } while (mA + mB > eCM);
  double s1    = m1 * m1;
  double s2    = m2 * m2;
  double sA    = mA * mA;
  double sB    = mB * mB; 

  // Energies and longitudinal momenta of excited hadrons.
  double eA    = 0.5 * (sCM + sA - sB) / eCM;
  double pzA   = sqrtpos(eA * eA - sA);
  Vec4   pA( 0., 0.,  pzA,       eA);
  Vec4   pB( 0., 0., -pzA, eCM - eA); 

  // Internal kinematics on side A, boost to CM frame and store constituents.
  if (excite1) { 
    double ec1   = 0.5 * (sA + mTsc1 - mTsac1) / mA; 
    double pzc1  = sqrtpos(ec1 * ec1 - mTsc1);
    if (rndmPtr->flat() > 0.5) pzc1 = -pzc1;  
    Vec4 pc1(   px1,  py1,  pzc1,      ec1);   
    Vec4 pac1( -px1, -py1, -pzc1, mA - ec1);   
    pc1.bst(pA);
    pac1.bst(pA);
    int newc = event.nextColTag();  
    event.append(  idc1, 63, i1, 0, 0, 0, newc,     0,  pc1,  mc1);
    event.append( idac1, 63, i1, 0, 0, 0, 0,     newc, pac1, mac1);
    event[i1].statusNeg();

  // Simple copy if not excited, ad set momentum as in collision frame.
  } else {
    int iNew = event.copy( i1, 63);
    event[iNew].p( pA);
    event[iNew].vProd( 0., 0., 0., 0.);
  }

  // Internal kinematics on side B, boost to CM frame and store constituents.
  if (excite2) { 
    double ec2   = 0.5 * (sB + mTsc2 - mTsac2) / mB; 
    double pzc2  = sqrtpos(ec2 * ec2 - mTsc2);
    if (rndmPtr->flat() > 0.5) pzc2 = -pzc2;  
    Vec4 pc2(   px2,  py2,  pzc2,      ec2);   
    Vec4 pac2( -px1, -py1, -pzc2, mA - ec2);   
    pc2.bst(pB);
    pac2.bst(pB);
    int newc = event.nextColTag();  
    event.append(  idc2, 63, i2, 0, 0, 0, newc,     0,  pc2,  mc2);
    event.append( idac2, 63, i2, 0, 0, 0, 0,     newc, pac2, mac2);
    event[i2].statusNeg();

  // Simple copy if not excited, ad set momentum as in collision frame.
  } else {
   int iNew = event.copy( i2, 63);
    event[iNew].p( pB);
    event[iNew].vProd( 0., 0., 0., 0.);
  }

  // Select t value and rotate outgoing particles accordingly.
  double lambda12 = pow2( sCM - s1 - s2) - 4. * s1 * s2;
  double lambdaAB = pow2( sCM - sA - sB) - 4. * sA * sB;
  double tLow     = -0.5 * (sCM - (s1 + s2 + sA + sB) + (s1 - s2) 
    * (sA - sB) / sCM + sqrtpos(lambda12 *  lambdaAB) / sCM);
  double tUpp     = ( (sA - s1) * (sB - s2) + (s1 + sB - s2 - sA)
    * (s1 * sB - s2 * sA) / sCM ) / tLow;
  double bNow     = bSlope();
  double eBtLow   = exp( bNow * tLow);
  double eBtUpp   = exp( bNow * tUpp);
  double tNow     = log( eBtLow + rndmPtr->flat() * (eBtUpp - eBtLow) ) / bNow; 
  double theta    = acos( (2. * tNow - tLow - tUpp) / (tUpp - tLow) );
  double phi      = 2. * M_PI * rndmPtr->flat();
  for (int i = sizeOld; i < event.size(); ++i) event[i].rot( theta, phi); 

  // Done.
  return true;
    
}

//-------------------------------------------------------------------------

// Split up hadron A into a colour-anticolour pair, with masses and pT values.
  
bool LowEnergyHadHad::splitA( double redMpT) {

  // Split up flavour of hadron into a colour and an anticolour constituent. 
  pair< int, int>  paircac  = splitFlav( id1 );
  idc1   = paircac.first;
  idac1  = paircac.second;
  if (idc1 == 0 || idac1 == 0) return false;

  // Find constituent masses and scale down to less than full mass.
  mc1    = particleDataPtr->m0( idc1);
  mac1   = particleDataPtr->m0( idac1);
  double redNow = redMpT * min( 1., m1 / (mc1 + mac1));
  mc1   *= redNow;
  mac1  *= redNow;

  // Select Gaussian relative transverse momenta for constituents.
  pair<double, double> gauss2 = rndmPtr->gauss2();
  px1    = redMpT * sigmaQ * gauss2.first;
  py1    = redMpT * sigmaQ * gauss2.second; 
  pTs1   = px1 * px1 + py1 * py1;

  // Construct transverse masses.
  mTsc1  = pow2(mc1)  + pTs1;
  mTsac1 = pow2(mac1) + pTs1;
  mTc1   = sqrt(mTsc1);
  mTac1  = sqrt(mTsac1); 

  // Done.
  return true;
    
}

//-------------------------------------------------------------------------

// Split up hadron B into a colour-anticolour pair, with masses and pT values.
  
bool LowEnergyHadHad::splitB( double redMpT) {

  // Split up flavour of hadron into a colour and an anticolour constituent. 
  pair< int, int>  paircac  = splitFlav( id2 );
  idc2   = paircac.first;
  idac2  = paircac.second;
  if (idc2 == 0 || idac2 == 0) return false;

  // Find constituent masses and scale down to less than full mass.
  mc2    = particleDataPtr->m0( idc2);
  mac2   = particleDataPtr->m0( idac2);
  double redNow = redMpT * min( 1., m2 / (mc2 + mac2));
  mc2   *= redNow;
  mac2  *= redNow;

  // Select Gaussian relative transverse momenta for constituents.
  pair<double, double> gauss2 = rndmPtr->gauss2();
  px2    = redMpT * sigmaQ * gauss2.first;
  py2    = redMpT * sigmaQ * gauss2.second; 
  pTs2   = px2 * px2 + py2 * py2;

  // Construct transverse masses.
  mTsc2  = pow2(mc2)  + pTs2;
  mTsac2 = pow2(mac2) + pTs2;
  mTc2   = sqrt(mTsc2);
  mTac2  = sqrt(mTsac2); 

  // Done.
  return true;
    
}

//-------------------------------------------------------------------------

// Split up a hadron into a colour and an anticolour part, of q or qq kinds.

pair< int, int> LowEnergyHadHad::splitFlav( int id) {
   
  // Hadron flavour content.
  int idAbs = abs(id);
  int iq1   = (idAbs/1000)%10;
  int iq2   = (idAbs/100)%10;
  int iq3   = (idAbs/10)%10;
  int iq4, iq5;

  // Nondiagonal mesons.
  if (iq1 == 0 && iq2 != iq3) {
    if (id != 130 && id != 310) {
      if (iq2%2 == 1) swap( iq2, iq3);
      if (id > 0) return make_pair( iq2, -iq3);
      else        return make_pair( iq3, -iq2);

    // K0S and K0L are mixes d sbar and dbar s.   
    } else {
      if (rndmPtr->flat() < 0.5) return make_pair( 3, -1);
      else                       return make_pair( 1, -3);
    }
      
  // Diagonal mesons: assume complete mixing ddbar and uubar.
  } else if (iq1 == 0) {
   if (iq2 < 3 || id == 331) {
     iq4 = (rndmPtr->flat() < 0.5) ? 1 : 2; 
     // eta and eta' can also be s sbar.
     if (id == 221 && rndmPtr->flat() < fracEtass) iq4 = 3;
     if (id == 331 && rndmPtr->flat() < fracEtaPss) iq4 = 3;
     return make_pair( iq4, -iq4);
   }

  // Octet baryons.
  } else if (idAbs%10 == 2) {
    // Three identical quarks: emergency in case of higher spin 1/2 multiplet.
    if (iq1 == iq2 && iq2 == iq3) {iq4 = iq1; iq5 = 1100 * iq1 + 3;}
    // Two identical quarks, like normal p or n.
    else if (iq1 == iq2 || iq2 == iq3) {
      double rr6 = 6. * rndmPtr->flat();
      if    (iq1 == iq2 && rr6 < 2.) { iq4 = iq3; iq5 = 1100 * iq1 + 3;}
      else if             (rr6 < 2.) { iq4 = iq1; iq5 = 1100 * iq3 + 3;}
      else if (rr6 < 3.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;}
      else               { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 1;}
    // Three nonidentical quarks, Sigma- or Lambda-like.
    } else { 
      int isp = (iq2 > iq3) ? 3 : 1;    
      if (iq3 > iq2) swap( iq2, iq3);   
      double rr12 = 12. * rndmPtr->flat();
      if      (rr12 < 4.) { iq4 = iq1; iq5 = 1000 * iq2 + 100 * iq3 + isp;}   
      else if (rr12 < 5.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + isp;}
      else if (rr12 < 6.) { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + isp;}
      else if (rr12 < 9.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 4 - isp;}
      else                { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + 4 - isp;}
    }
    if (id > 0) return make_pair(  iq4,  iq5);
    else        return make_pair( -iq5, -iq4);
    
  // Decuplet baryons.
  } else {
    double rr3 = 3. * rndmPtr->flat();
    if (rr3 < 1.)      { iq4 = iq1; iq5 = 1000 * iq2 + 100 * iq3 + 3;}
    else if (rr3 < 2.) { iq4 = iq2; iq5 = 1000 * iq1 + 100 * iq3 + 3;} 
    else               { iq4 = iq3; iq5 = 1000 * iq1 + 100 * iq2 + 3;} 
    if (id > 0) return make_pair(  iq4,  iq5);
    else        return make_pair( -iq5, -iq4); 
  }

  // Done. (Fake call to avoid unwarranted compiler warning.)
  return make_pair( 0, 0);

}

//-------------------------------------------------------------------------

// Find relative momentum of colour and anticolour constituents in hadron.

double LowEnergyHadHad::splitZ( int iq1, int iq2, double mRat1, double mRat2) {

  // Initial values.
  int iq1Abs = abs(iq1);
  int iq2Abs = abs(iq2);
  if (iq2Abs > 10) swap( mRat1, mRat2);
  double x1, x2, x1a, x1b;

  // Handle mesons.
  if (iq1Abs < 10 && iq1Abs < 10) {
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
    if (iq2Abs > 10) swap( x1, x2);
  }

  // Return z value.      
  return x1 / (x1 + x2); 
  
}

//-------------------------------------------------------------------------

// Pick slope b of exp(b * t) for elastic and diffractive events.

double LowEnergyHadHad::bSlope() {
  return 5.;
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
  pythia.readString("Next:numberShowEvent = 0");
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

   // Find an original hadron pair with invariant mass above 10 GeV.
   int ncollide = 0;
   int sizeOrig = event.size();
   cout << " original event size is " << sizeOrig << endl;
   for (int i1 = 0; i1 < sizeOrig - 1; ++i1) 
   if (event[i1].isFinal() && event[i1].isHadron()) {
     for (int i2 = i1 + 1; i2 < sizeOrig; ++i2) 
     if (event[i2].isFinal() && event[i2].isHadron() 
     && (event[i1].p() + event[i2].p()).mCalc() > 10.) {
       //int kind = 1. + 5. * pythia.rndm.flat();
       int kind = 3;
       cout << " collision between " << i1 << " and " << i2 
            << " is of kind " << kind << endl;
       LEHH.collide( kind, i1, i2, event);
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
