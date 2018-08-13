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

// RescatteringEvent data structure

//--------------------------------------------------------------------------

class RescatteringEvent {
public:
// @TODO
  RescatteringEvent(int iDecayIn, Vec4 originIn)
    : iFirst(iDecayIn), iSecond(0), origin(originIn) {}
  RescatteringEvent(int iFirstIn, int iSecondIn, Vec4 originIn)
    : iFirst(iFirstIn), iSecond(iSecondIn), origin(originIn) {}
  
  bool isDecay() { return iSecond == 0; }

  int iFirst, iSecond;
  Vec4 origin;
};

//--------------------------------------------------------------------------

struct RescatteringEventComparer {
  bool operator()(const RescatteringEvent& lhs, const RescatteringEvent& rhs) {
    return lhs.origin.e() > rhs.origin.e();
  }
};

//==========================================================================

// The Rescattering class

//--------------------------------------------------------------------------

bool Rescattering::calculateDecay(Particle& pIn, Vec4& originOut) {
  // @TODO: Also check maximum lifetime
  if (!pIn.canDecay() || !pIn.mayDecay() 
      || pIn.tau() > tau0Max)
    return false;

  originOut = pIn.vDec();
  return true;
}

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

bool Rescattering::calculateInteraction(int idA, int idB, 
  Event& event, Vec4& originOut) {

  // Abort if the two particles come from the same decay
  // @TODO: Test what happens if this gets turned off
  if (event[idA].mother1() == event[idB].mother1()
      && (event[idA].status() >= 90 && event[idA].status() <= 99))
    return false;

  // Get cross section from data
  // @TODO: double sigma = crossSectionDataPtr->sigma(hadA.id(), hadB.id());
  // if (sigma == 0) return false;
  double sigma = 40;

  // Set up positions for each particle in their CM frame
  RotBstMatrix frame;
  frame.toCMframe(event[idA].p(), event[idB].p());

  Vec4 vA = event[idA].vProd();
  Vec4 pA = event[idA].p();
  Vec4 vB = event[idB].vProd();
  Vec4 pB = event[idB].p();

  vA.rotbst(frame); vB.rotbst(frame);
  pA.rotbst(frame); pB.rotbst(frame);

  // Abort if impact parameter is too large
  if ((vA - vB).pT2() > (0.1 * FM2MM * FM2MM) * sigma / M_PI)
    return false;

  // Check if particles have already passed each other
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

  // @TODO: Abort if final separation is spacelike?

  // Return event candidate
  originOut = origin;
  return true;
}

//--------------------------------------------------------------------------

static void phaseSpace2(Rndm* rndm, Vec4 pTotIn, double mA, double mB,
                        Vec4& p1Out, Vec4& p2Out) {
  
  double m0 = pTotIn.mCalc();

  double eA = 0.5 * (m0 + (mA * mA - mB * mB) / m0);
  double eB = 0.5 * (m0 + (mB * mB - mA * mA) / m0);

  double p = 0.5 / m0 * sqrtpos((m0 + (mA + mB)) * (m0 - (mA + mB))
                              * (m0 + (mA - mB)) * (m0 - (mA - mB))); 
															
  double phi = 2 * M_PI * rndm->flat();
  double costh = 2 * rndm->flat() - 1;
  double sinth = sqrt(1 - costh * costh);

  double x = p * sinth * cos(phi),
         y = p * sinth * sin(phi),
         z = p * costh;

  p1Out = Vec4(x, y, z, eA);
  p1Out.bst(pTotIn);
  p2Out = Vec4(-x, -y, -z, eB);
  p2Out.bst(pTotIn);

  double err = (p1Out + p2Out - pTotIn).pAbs() 
                / (p1Out + p2Out + pTotIn).pAbs();
  if (isnan(err) || isinf(err) || err > 1.0e-9)
  {
    cout << "Phase space not quite good " << endl
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

  int status = (hadA.status() == 111 || hadB.status() == 111) ? 112 : 111;

  Vec4 mom1, mom2;
  phaseSpace2(this->rndmPtr, hadA.p() + hadB.p(),
              hadA.m0(), hadB.m0(),
              mom1, mom2);
  
  vector<Particle> newParticles;

  newParticles.push_back(Particle(hadA.id(), status, idA, idB, 0, 0, 
                                  0, 0, mom1, mom1.mCalc()));
  newParticles.push_back(Particle(hadB.id(), status, idA, idB, 0, 0, 
                                  0, 0, mom2, mom2.mCalc()));

  // @TODO Some value copying going on here, but this is placeholder anyway
  for (auto particle : newParticles)
    event.append(particle);

  // Update the interacting particles
  for (int i : { idA, idB }) {

    event[i].statusNeg();
    event[i].daughters(oldSize, event.size() - 1);

    // Set proper lifetime (decay vertex the point closest to origin)
    // @TODO: Test that the lifetime is set correctly
    Vec4 beta = event[i].p() / event[i].e();
    double gamma = 1 / sqrt(1 - dot3(beta, beta));
    double t = dot3(origin - event[i].vProd(), beta) / dot3(beta, beta);
    event[i].tau(t / gamma);
  }
}

//-------------------------------------------------------------------------- 

bool Rescattering::canScatter(Particle& particle)
{
  return particle.isFinal() && particle.isHadron()
      && particle.vProd().pAbs() < radiusMax;
}

void Rescattering::next(Event& event) {
  
  auto candidates = std::priority_queue<RescatteringEvent,
                                        vector<RescatteringEvent>,
                                        RescatteringEventComparer>();

  for (int iFirst = 0; iFirst < event.size(); ++iFirst) {
    if (!canScatter(event[iFirst]))
      continue;

    Vec4 origin;
    if (calculateDecay(event[iFirst], origin)) {
      candidates.push(RescatteringEvent(iFirst, origin));
    }

    for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
      if (!canScatter(event[iSecond]))
        continue; 
      if (calculateInteraction(iFirst, iSecond, event, origin)) {
        candidates.push(RescatteringEvent(iFirst, iSecond, origin));
      }
    }
  }

  while (!candidates.empty()) {
    RescatteringEvent ev = candidates.top();
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
    {
      for (int iFirst = oldSize; iFirst < event.size(); ++iFirst) {
        if (!canScatter(event[iFirst]))
          continue;
        
        Vec4 origin;
        if (calculateDecay(event[iFirst], origin))
          candidates.push(RescatteringEvent(iFirst, origin));

        for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
          if (!canScatter(event[iFirst]))
            continue;

          if (calculateInteraction(iFirst, iSecond, event, origin))
            candidates.push(RescatteringEvent(iFirst, iSecond, origin));
        }
      }
    }
  }
}

} // end namespace Pythia8
