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
			|| pIn.tau() > settingsPtr->parm("Rescattering:tau0Max"))
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

bool Rescattering::calculateInteraction(Particle& p1In,
	Particle& p2In, Vec4& originOut) {

	// Order particles so that p1 is produced before p2
	Particle& p1 = p1In.tProd() < p2In.tProd() ? p1In : p2In;
	Particle& p2 = p1In.tProd() < p2In.tProd() ? p2In : p1In;

	if ((p1.statusAbs() == 111 || p2.statusAbs() == 111)
		  && !settingsPtr->flag("Rescattering:allowSecondRescattering"))
		return false;

	// Get cross section from data
  double sigma = crossSectionDataPtr->sigma(p1.id(), p2.id());

	// Abort if the particles cannot interact
	if (sigma == 0)
		return false;

	// Set up positions and velocities for each particle in their CM frame
	RotBstMatrix m;
	m.toCMframe(p1.p(), p2.p());

	Vec4 v1 = p1.p() / p1.e(); 
	Vec4 v2 = p2.p() / p2.e();
	Vec4 x1 = p1.vProd(); 
	Vec4 x2 = p2.vProd();
	v1.rotbst(m); v2.rotbst(m);
	x1.rotbst(m); x2.rotbst(m);
		
	// Move the first particle to time origin of the second particle
	x1 += v1 * (x2.e() - x1.e());

	// Abort if the particles are not moving towards each other
	if (dot3(x2 - x1, v2 - v1) >= 0)
		return false;

	// Abort if the impact parameter is too large
	if ((x2 - x1).pT() > sqrt(sigma / M_PI))
		return false;

	// Calculate origin in CM frame and transform to lab frame
	double dv = v1.pz() - v2.pz();
	double tColl = -(x1 - x2).pz() / dv;
	Vec4 collisionOrigin((x1.px() + x2.px()) / 2, (x1.py() + x2.py()) / 2,
										   x1.pz() + tColl * v1.pz(), x1.e() + tColl);
	m.invert();
	collisionOrigin.rotbst(m);

	if ((userHooksPtr != nullptr) && userHooksPtr->canVetoRescatteringInteraction())
		if (userHooksPtr->doVetoRescatteringInteraction(p1, p2, collisionOrigin))
			return false;

	// Return event candidate
	originOut = collisionOrigin;
	return true;
}

//--------------------------------------------------------------------------

static void phaseSpace2(Rndm* rndm, Vec4 pTotIn, double mAIn, double mBIn,
											  Vec4& p1Out, Vec4& p2Out) {
	double phi = 2 * M_PI * rndm->flat();
	double theta = acos(2 * rndm->flat() - 1);

	double mA2 = mAIn * mAIn, mB2 = mBIn * mBIn;
	double s = pTotIn.m2Calc();
	double p = sqrt(pow2(mA2 - mB2) / (4. * s) 
									- (mA2 + mB2) / 2. + s / 4.);
	
	double x = p * cos(phi) * sin(theta),
				 y = p * sin(phi) * sin(theta),
				 z = p * cos(theta);

	p1Out = Vec4(x, y, z, sqrt(p * p + mA2));
	p1Out.bst(pTotIn);
	p2Out = Vec4(-x, -y, -z, sqrt(p * p + mB2));
	p2Out.bst(pTotIn);
}

void Rescattering::produceScatteringProducts(int id1, int id2, 
	Vec4& origin, Event& event) {

	Particle& p1 = event[id1];
	Particle& p2 = event[id2]; 
	int oldSize = event.size();
	
  // Pick the interaction channel
	CrossSectionDataEntry* entry 
		= crossSectionDataPtr->findCrossSection(p1.id(), p2.id());
	InteractionChannel channel = entry->pickChannel();

	int status = (p1.status() == 111 || p2.status() == 111) ? 112 : 111;

	Vec4 mom1, mom2;
	phaseSpace2(this->rndmPtr, p1.p() + p2.p(),
							particleDataPtr->m0(channel.product(0)),
							particleDataPtr->m0(channel.product(1)),
				      mom1, mom2);
	
	// @TODO: More than two particles
	event.append(channel.product(0), status, id1, id2, 0, 0, 0, 0,
	  					 mom1, mom1.mCalc());
	event.append(channel.product(1), status, id1, id2, 0, 0, 0, 0,
	  					 mom2, mom2.mCalc());

	// Update the interacting particles
	for (int i : { id1, id2 }) {

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

void Rescattering::next(Event& event) {

	auto candidates = std::priority_queue<RescatteringEvent,
																				vector<RescatteringEvent>,
																				RescatteringEventComparer>();
	
	for (int iFirst = 0; iFirst < event.size(); ++iFirst) {
		if (!event[iFirst].isFinal())
			continue;

		Vec4 origin;
		if (calculateDecay(event[iFirst], origin)) {
			candidates.push(RescatteringEvent(iFirst, origin));
		}

		for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
			if (!event[iSecond].isFinal())
				continue; 
			if (calculateInteraction(event[iFirst], event[iSecond], origin)) {
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
		for (int iFirst = oldSize; iFirst < event.size(); ++iFirst) {
			if (!event[iFirst].isFinal())
				continue;
			
			Vec4 origin;
			if (calculateDecay(event[iFirst], origin))
				candidates.push(RescatteringEvent(iFirst, origin));

			for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
				if (!event[iSecond].isFinal())
					continue;

				if (calculateInteraction(event[iFirst], event[iSecond], origin))
					candidates.push(RescatteringEvent(iFirst, iSecond, origin));
			}
		}
	}
}

} // end namespace Pythia8
