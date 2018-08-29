#include <algorithm>
#include <functional>
#include <initializer_list>
#include <queue>
#include "Pythia8/CrossSectionData.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/Rescattering.h"

 
//--------------------------------------------------------------------------

namespace Pythia8 {

//==========================================================================

// PriorityVertex data structure

//--------------------------------------------------------------------------

struct Rescattering::PriorityVertex {
public:
// @TODO
  PriorityVertex(int iDecayIn, Vec4 originIn)
    : iFirst(iDecayIn), iSecond(0), origin(originIn) {}
  PriorityVertex(int iFirstIn, int iSecondIn, Vec4 originIn)
    : iFirst(iFirstIn), iSecond(iSecondIn), origin(originIn) {}
  
  bool isDecay() { return iSecond == 0; }

  // Priority comparison; the one with lowest time has highest priority!
  bool operator<(const PriorityVertex& r) const
  { return origin.e() > r.origin.e(); }

  int iFirst, iSecond;
  Vec4 origin;
};

//==========================================================================

// The Rescattering class

//--------------------------------------------------------------------------

bool Rescattering::produceDecayProducts(int iDec, Event& event) {

  int oldSize = event.size();

  // @TODO: Something with the return value
  decays.decay(iDec, event);

  for (int i = oldSize; i < event.size(); ++i)
    event[i].status(113);

  return true;
}



//--------------------------------------------------------------------------

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
  if (isnan(err) || isinf(err) || err > 1.0e-9)
  {
    cerr << "Phase space not quite good " << endl
         << "   In: " << pTotIn 
         << "  Out: " << p1Out + p2Out
         << "error: " << scientific << err << endl
         << endl;
  }

}

void Rescattering::produceScatteringProducts(int idA, int idB, 
  Vec4& origin, Event& event) {

  Particle& hadA = event[idA];
  Particle& hadB = event[idB]; 
  int oldSize = event.size();

  /* @TODO
  // Pick the interaction channel
  CrossSectionDataEntry* entry 
    = crossSectionDataPtr->findCrossSection(hadA.id(), hadB.id());

  if (entry == nullptr) 
  {
    // @TODO: Output error
    return;
  }
  
  InteractionChannel channel = entry->pickChannel();

  int status = (hadA.status() == 111 || hadB.status() == 111) ? 112 : 111;

  Vec4 mom1, mom2;
  phaseSpace2(this->rndmPtr, hadA.p() + hadB.p(),
              particleDataPtr->m0(channel.product(0)),
              particleDataPtr->m0(channel.product(1)),
              mom1, mom2);
  
  vector<Particle> newParticles;

  newParticles.push_back(Particle(channel.product(0), status, idA, idB, 0, 0, 
                                  0, 0, mom1, mom1.mCalc()));
  newParticles.push_back(Particle(channel.product(1), status, idA, idB, 0, 0, 
                                  0, 0, mom2, mom2.mCalc()));
*/

  int status = (hadA.status() == 111 || hadA.status() == 112
             || hadB.status() == 111 || hadB.status() == 112) ? 112 : 111;

  Vec4 mom1, mom2;
  phaseSpace2(this->rndmPtr, hadA.p() + hadB.p(),
              hadA.m(), hadB.m(),
              mom1, mom2);
  
  vector<Particle> newParticles;

  newParticles.push_back(Particle(hadA.id(), status, idA, idB, 0, 0, 
                                  0, 0, mom1, mom1.mCalc()));
  newParticles.push_back(Particle(hadB.id(), status, idA, idB, 0, 0, 
                                  0, 0, mom2, mom2.mCalc()));

  // @TODO Some value copying going on here, but this is placeholder anyway
  for (auto particle : newParticles)
  {
    particle.vProd(origin);
    event.append(particle);
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
}

//-------------------------------------------------------------------------- 

void Rescattering::calcDecaysRescatters(Event& event, int iStart,
                            priority_queue<PriorityVertex>& queue)
{
  // @TODO Stress test and optimise if needed
  for (int iFirst = iStart; iFirst < event.size(); ++iFirst) {
    Particle& hadA = event[iFirst];

    if (!hadA.isFinal() || !hadA.isHadron() 
        || hadA.vProd().pAbs() >= radiusMax)
      continue;

    if (doDecays && hadA.canDecay() && hadA.mayDecay())
      queue.push(PriorityVertex(iFirst, hadA.vDec()));

    for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
      Particle& hadB = event[iSecond];
      if (!hadB.isFinal() || !hadB.isHadron() 
          || hadB.vProd().pAbs() >= radiusMax)
        continue;
      
      // Continue if the two particles come from the same decay
      // @TODO: Test what happens if this gets turned off
      if (hadA.mother1() == hadB.mother1() 
          && (hadA.status() >= 90 && hadA.status() <= 99))
        continue;

      // Get cross section from data
      // @TODO: actually get cross section from somewhere
      double sigma = 40;

      // @TODO: Ideally, we just care about the invariant closest distance and
      //  the time of closest approach at this point. All these calculations
      //  could be shortened. In particular, profiling shows that
      //  frame.toCMframe takes a significant part of the running time

      // Set up positions for each particle in their CM frame
      RotBstMatrix frame;
      frame.toCMframe(hadA.p(), hadB.p());

      Vec4 vA = hadA.vProd();
      Vec4 pA = hadA.p();
      Vec4 vB = hadB.vProd();
      Vec4 pB = hadB.p();

      vA.rotbst(frame); vB.rotbst(frame);
      pA.rotbst(frame); pB.rotbst(frame);

      // Abort if impact parameter is too large
      if ((vA - vB).pT2() > MB2MMSQ * sigma / M_PI)
        continue;

      // Check if particles have already passed each other
      double t0 = max(vA.e(), vB.e());
      double zA = vA.pz() + (t0 - vA.e()) * pA.pz() / pA.e();
      double zB = vB.pz() + (t0 - vB.e()) * pB.pz() / pB.e();

      if (zA >= zB)
        continue;

      // Calculate collision origin and transform to lab frame
      double tCollision = t0 - (zB - zA) / (pB.pz() / pB.e() - pA.pz() / pA.e());
      Vec4 origin(0.5 * (vA.px() + vB.px()),
                  0.5 * (vA.py() + vB.py()),
                  zA + pA.pz() / pA.e() * (tCollision - t0),
                  tCollision);

      frame.invert();
      origin.rotbst(frame);

      // Return event candidate
      queue.push(PriorityVertex(iFirst, iSecond, origin));
    }
  }
}

void Rescattering::next(Event& event) {
  
  priority_queue<PriorityVertex> candidates; 

  calcDecaysRescatters(event, 0, candidates);

  while (!candidates.empty()) {
    PriorityVertex ev = candidates.top();
    candidates.pop();

    // Abort if either particle has already interacted elsewhere
    if (!event[ev.iFirst].isFinal() 
    || (!ev.isDecay() && !event[ev.iSecond].isFinal()))
      continue;

    int oldSize = event.size();

    // Produce products
    if (ev.isDecay()) 
      produceDecayProducts(ev.iFirst, event);
    else
      produceScatteringProducts(ev.iFirst, ev.iSecond, ev.origin, event);

    // Check for new interactions
    if (doSecondRescattering)
      calcDecaysRescatters(event, oldSize, candidates);
  }
}

} // end namespace Pythia8
