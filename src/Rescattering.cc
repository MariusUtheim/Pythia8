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

// The Rescattering class

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

void Rescattering::rescatter(int idA, int idB, 
  Vec4 origin, Event& event) {

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

} // end namespace Pythia8
