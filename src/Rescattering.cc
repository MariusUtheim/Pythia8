#include <algorithm>
#include <functional>
#include <initializer_list>
#include <queue>

#include "Pythia8/Rescattering.h"

 
//--------------------------------------------------------------------------

namespace Pythia8 {

//==========================================================================

// The Rescattering class

//--------------------------------------------------------------------------

void Rescattering::init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn, ParticleData* particleDataPtrIn)
  {
    infoPtr = infoPtrIn;
    rndmPtr = rndmPtrIn; 
    particleDataPtr = particleDataPtrIn;
    
    leHadHad.init(infoPtrIn, settings, particleDataPtrIn, rndmPtrIn);
}

// Temporary function for sampling 2D phase space
static void phaseSpace2(Rndm* rndmPtr, Vec4 pTotIn, double mA, double mB,
                        Vec4& p1Out, Vec4& p2Out) {
  
  double m0 = pTotIn.mCalc();

  double eA = 0.5 * (m0 + (mA * mA - mB * mB) / m0);
  double eB = 0.5 * (m0 + (mB * mB - mA * mA) / m0);

  double p = (0.5 / m0) * sqrtpos((m0 + (mA + mB)) * (m0 - (mA + mB))
                                * (m0 + (mA - mB)) * (m0 - (mA - mB))); 

  double phi = 2. * M_PI * rndmPtr->flat();
  double costh = 2. * rndmPtr->flat() - 1.;
  double sinth = sqrt(1 - costh * costh);

  double px = p * sinth * cos(phi),
         py = p * sinth * sin(phi),
         pz = p * costh;

  p1Out = Vec4(px, py, pz, eA);
  p1Out.bst(pTotIn);
  p2Out = Vec4(-px, -py, -pz, eB);
  p2Out.bst(pTotIn);
 
  double err = (p1Out + p2Out - pTotIn).pAbs() 
                / (p1Out + p2Out + pTotIn).pAbs();
  if (isnan(err) || isinf(err) || err > 1.0e-8)
  {
    cerr << "Phase space not quite good " << endl
         << "   In: " << pTotIn 
         << "  Out: " << p1Out + p2Out
         << "error: " << scientific << err << endl
         << endl;
  }

}

static vector<Vec4> phaseSpace(Rndm* rndmPtr, Vec4 pTotIn, vector<double> msIn) {
  Vec4 p1, p2;
  phaseSpace2(rndmPtr, pTotIn, msIn[0], msIn[1], p1, p2);
  return vector<Vec4> { p1, p2 };
}

void Rescattering::rescatter(int idA, int idB, 
  Vec4 origin, Event& event) {

  leHadHad.collide(idA, idB, 0, event, origin);

/*
  Particle& hadA = event[idA];
  Particle& hadB = event[idB]; 

  vector<int> products = { hadA.id(), hadB.id() };//= crossSecPtr->pickProducts(hadA.id(), hadB.id(), eCM);
 
  int oldSize = event.size();
  int status = (hadA.status() == 111 || hadA.status() == 112
             || hadB.status() == 111 || hadB.status() == 112) ? 112 : 111;

  vector<double> masses(products.size());
  for (int i = 0; i < products.size(); ++i)
    masses[i] = particleDataPtr->m0(products[i]);

  vector<Vec4> momenta = phaseSpace(rndmPtr, hadA.p() + hadB.p(), masses);

  for (int i = 0; i < products.size(); ++i)
  {
    Particle newParticle(products[i], status, idA, idB, 0, 0,
      0, 0, momenta[i], momenta[i].mCalc());
    newParticle.vProd(origin);
    event.append(newParticle);
  }

  // Update the interacting particles
  for (int i : { idA, idB }) {

    event[i].statusNeg();
    event[i].daughters(oldSize, event.size() - 1);

    // Set proper lifetime (decay vertex the point closest to origin)
    // @TODO: Test that the lifetime is set correctly
    Vec4 beta = event[i].p() / event[i].e();
    double gamma = event[i].e() / event[i].m();
    double t = dot3(origin - event[i].vProd(), beta) / dot3(beta, beta);
    event[i].tau(t / gamma);
  }
  */
}

//--------------------------------------------------------------------------

bool Rescattering::calcRescatterOrigin(int idA, int idB, Event& event, 
  Vec4& originOut)
{
  Particle& hadA = event[idA];
  Particle& hadB = event[idB];

  // @TODO: Profiling shows that frame.toCMframe is the most significant
  // bottleneck in Pythia for high-multiplicity events. We should think about
  // checks that can be made to abort early (e.g. particles moving away
  // from each other in the lab frame)

  // Set up positions for each particle in their CM frame
  RotBstMatrix frame;
  frame.toCMframe(hadA.p(), hadB.p());

  Vec4 vA = hadA.vProd();
  Vec4 pA = hadA.p();
  Vec4 vB = hadB.vProd();
  Vec4 pB = hadB.p();

  vA.rotbst(frame); vB.rotbst(frame);
  pA.rotbst(frame); pB.rotbst(frame);

  double eCM = (pA + pB).mCalc();
  double sigma = leHadHad.sigmaTotal(hadA.idAbs(), hadB.idAbs(), eCM);

  // Abort if impact parameter is too large
  if ((vA - vB).pT2() > MB2MMSQ * sigma / M_PI)
    return false;

  // Abort if particles have already passed each other
  double t0 = max(vA.e(), vB.e());
  double zA = vA.pz() + (t0 - vA.e()) * pA.pz() / pA.e();
  double zB = vB.pz() + (t0 - vB.e()) * pB.pz() / pB.e();

  if (zA >= zB)
    return false;

  // Calculate collision origin and transform to lab frame
  double tCollision = t0 - (zB - zA) / (pB.pz() / pB.e() - pA.pz() / pA.e());
  Vec4 origin(0.5 * (vA.px() + vB.px()),
              0.5 * (vA.py() + vB.py()),
              zA + pA.pz() / pA.e() * (tCollision - t0),
              tCollision);

  frame.invert();
  origin.rotbst(frame);

  // Return 
  originOut = origin;
  return true;
}

} // end namespace Pythia8
