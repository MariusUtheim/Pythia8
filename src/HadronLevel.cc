// HadronLevel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the HadronLevel class.

#include "Pythia8/HadronLevel.h"

namespace Pythia8 {

//==========================================================================

// The HadronLevel class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.

// Small safety mass used in string-end rapidity calculations.
const double HadronLevel::MTINY = 0.1;

//--------------------------------------------------------------------------

// Node for ordering scatterings and decays

struct HadronLevel::PriorityNode {
// @TODO
  PriorityNode(int iDecayIn, Vec4 originIn)
    : iFirst(iDecayIn), iSecond(0), origin(originIn) {}
  PriorityNode(int iFirstIn, int iSecondIn, Vec4 originIn)
    : iFirst(iFirstIn), iSecond(iSecondIn), origin(originIn) {}
  
  bool isDecay() { return iSecond == 0; }

  // Priority comparison to be used by priority_queue
  // Note that lower t means higher priority!
  bool operator<(const PriorityNode& r) const
  { return origin.e() > r.origin.e(); }

  int iFirst, iSecond;
  Vec4 origin;
};

//--------------------------------------------------------------------------

// Find settings. Initialize HadronLevel classes as required.

bool HadronLevel::init(Info* infoPtrIn, Settings& settings,
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn,
  Couplings* couplingsPtrIn, TimeShower* timesDecPtr,
  RHadrons* rHadronsPtrIn, DecayHandler* decayHandlePtr,
  vector<int> handledParticles, UserHooks* userHooksPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  particleDataPtr = particleDataPtrIn;
  rndmPtr         = rndmPtrIn;
  couplingsPtr    = couplingsPtrIn;
  rHadronsPtr     = rHadronsPtrIn;
  userHooksPtr    = userHooksPtrIn;

  // Main flags.
  doHadronize     = settings.flag("HadronLevel:Hadronize");
  doDecay         = settings.flag("HadronLevel:Decay");
  doRescatter     = settings.flag("HadronLevel:Rescatter");
  doBoseEinstein  = settings.flag("HadronLevel:BoseEinstein");

  // Boundary mass between string and ministring handling.
  mStringMin      = settings.parm("HadronLevel:mStringMin");

  // For junction processing.
  eNormJunction   = settings.parm("StringFragmentation:eNormJunction");

  // Allow R-hadron formation.
  allowRH         = settings.flag("RHadrons:allow");

  // Allow decayed/rescattered particles to rescatter again
  allowMultipleRescatterings
                  = settings.flag("Rescattering:allowMultipleRescatterings");

  // Particles that should decay or not before Bose-Einstein stage.
  widthSepBE      = settings.parm("BoseEinstein:widthSep");

  // Need string density information be collected?
  closePacking    = settings.flag("StringPT:closePacking");

  // Rope hadronization. Setting of partonic production vertices.
  doRopes         = settings.flag("Ropewalk:RopeHadronization");
  doShoving       = settings.flag("Ropewalk:doShoving");
  doFlavour       = settings.flag("Ropewalk:doFlavour");
  doVertex        = settings.flag("PartonVertex:setVertex");
  doBuffon        = settings.flag("Ropewalk:doBuffon");

  // Initialize Ropewalk and Flavour Ropes.
  if (doRopes) {
    if (!ropewalk.init(infoPtr, settings, rndmPtr)) return false;
    flavourRope.init(&settings, rndmPtr, particleDataPtr, infoPtr, &ropewalk);
  }

  // Initialize auxiliary fragmentation classes.
  flavSel.init(settings,  particleDataPtr, rndmPtr, infoPtr);
  pTSel.init(  settings,  particleDataPtr, rndmPtr, infoPtr);
  zSel.init(   settings, *particleDataPtr, rndmPtr, infoPtr);

  // Initialize auxiliary administrative class.
  colConfig.init(infoPtr, settings, &flavSel);

  // Initialize string and ministring fragmentation.
  stringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr,
    &flavSel, &pTSel, &zSel, &flavourRope, userHooksPtr);
  ministringFrag.init(infoPtr, settings, particleDataPtr, rndmPtr,
    &flavSel, &pTSel, &zSel);

  // Initialize particle decays and rescatterings.
  decays.init(infoPtr, settings, particleDataPtr, rndmPtr, couplingsPtr,
    timesDecPtr, &flavSel, decayHandlePtr, handledParticles);

  rescatterings.init(infoPtr, particleDataPtr, rndmPtr,
    nullptr /* @TODO Replace by CrossSectionData */, userHooksPtr);

  if (doRescatter && !settings.flag("Fragmentation:setVertices"))
  {
    // @TODO Handle error message correctly and at the right location
    infoPtr->errorMsg("Error in HadronLevel::init: HadronLevel:Rescatter "
      "is on, but Fragmentation:setVertices is off.");
  }

  // Initialize BoseEinstein.
  boseEinstein.init(infoPtr, settings, *particleDataPtr);

  // Initialize Hidden-Valley fragmentation, if necessary.
  useHiddenValley = hiddenvalleyFrag.init(infoPtr, settings,
    particleDataPtr, rndmPtr);

  // Send flavour and z selection pointers to R-hadron machinery.
  rHadronsPtr->fragPtrs( &flavSel, &zSel);

  // Initialize the colour tracing class.
  colTrace.init(infoPtr);

  // Initialize the junction splitting class.
  junctionSplitting.init(infoPtr, settings, rndmPtr, particleDataPtr);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Hadronize and decay the next parton-level.

bool HadronLevel::next( Event& event) {

  // Store current event size to mark Parton Level content.
  event.savePartonLevelSize();

  // Do Hidden-Valley fragmentation, if necessary.
  if (useHiddenValley) hiddenvalleyFrag.fragment(event);

  // Colour-octet onia states must be decayed to singlet + gluon.
  if (!decayOctetOnia(event)) return false;

  // Set lifetimes for already existing hadrons, like onia.
  for (int i = 0; i < event.size(); ++i) if (event[i].isHadron())
    event[i].tau( event[i].tau0() * rndmPtr->exp() );

  // Remove junction structures.
  if (!junctionSplitting.checkColours(event)) {
    infoPtr->errorMsg("Error in HadronLevel::next: "
        "failed colour/junction check");
    return false;
  }

  // Possibility of hadronization inside decay, but then no BE second time.
  bool decaysProducedPartons;
  bool doBoseEinsteinNow = doBoseEinstein;
  do {
    decaysProducedPartons = false;

    // First part: string fragmentation.
    if (doHadronize) {

      // Find the complete colour singlet configuration of the event.
      // Keep junctions if we do shoving.
      if (!findSinglets( event, (doRopes && doShoving) )) return false;

      // Fragment off R-hadrons, if necessary.
      if (allowRH && !rHadronsPtr->produce( colConfig, event))
        return false;

      // Save list with rapidity pairs of the different string pieces.
      if (closePacking) {
        vector< vector< pair<double,double> > > rapPairs =
          rapidityPairs(event);
        colConfig.rapPairs = rapPairs;
      }

      // Let strings interact in rope hadronization treatment.
      if (doRopes) {

        // Do the shoving treatment.
        if (doShoving) {
          // For shoving we need explicit vertex information.
          if (!doVertex) {
            infoPtr->errorMsg("Error in HadronLevel::next: "
              "shoving enabled, but no vertex info.");
            return false;
          }
          // Extract all string segments from the event.
          ropewalk.extractDipoles(event, colConfig);
          // String shoving.
          ropewalk.shoveTheDipoles(event);
          // Find singlets again.
          iParton.resize(0);
          colConfig.clear();
          if (!findSinglets( event)) {
            infoPtr->errorMsg("Error in HadronLevel::next: "
              "ropes: failed 2nd singlet tracing.");
            return false;
          }
        }

        // Prepare for flavour ropes.
        if (doFlavour) {
          if (doVertex && !doBuffon) {
            ropewalk.extractDipoles(event, colConfig);
            ropewalk.calculateOverlaps();
          }
          // Else default to Buffon treatment which
          // does not need dipole extraction and overlaps.
          else {
            infoPtr->errorMsg("Error in HadronLevel::next: "
              "ropes: Flavour enabled, but no space time information.");
          }
        }
      }

      // Process all colour singlet (sub)systems.
      for (int iSub = 0; iSub < colConfig.size(); ++iSub) {

        // Collect sequentially all partons in a colour singlet subsystem.
        colConfig.collect(iSub, event);

        // String fragmentation of each colour singlet (sub)system.
        if ( colConfig[iSub].massExcess > mStringMin ) {
          if (!stringFrag.fragment( iSub, colConfig, event)) return false;

        // Low-mass string treated separately. Tell if diffractive system.
        } else {
          bool isDiff = infoPtr->isDiffractiveA() || infoPtr->isDiffractiveB();
          if (!ministringFrag.fragment( iSub, colConfig, event, isDiff))
            return false;
        }
      }
    }

    // Second part: sequential decays of short-lived particles (incl. K0).

    // If rescattering is off, we don't care about the order of the decays
    if (doDecay && !doRescatter) {
      decaysProducedPartons = decays.decayAll(event, widthSepBE);
    }
    // If rescattering is on, decays/rescatterings must happen in order
    else if (doRescatter) {
      
      //@TODO: Pick the best underlying data structure (probably red-black tree)
      priority_queue<PriorityNode> candidates; 

      queueDecaysAndRescatterings(event, 0, candidates);

      while (!candidates.empty()) {
        PriorityNode node = candidates.top();
        candidates.pop();

        // Abort if either particle has already interacted elsewhere
        if (!event[node.iFirst].isFinal() 
        || (!node.isDecay() && !event[node.iSecond].isFinal()))
          continue;

        int oldSize = event.size();

        // Perform the queued action
        if (node.isDecay()) {
          decays.decay(node.iFirst, event);
          if (decays.moreToDo()) decaysProducedPartons = true;
        }
        else
          rescatterings.rescatter(node.iFirst, node.iSecond, node.origin, event);

        // Check for new interactions
        // @TODO Allow decays of rescattered hadrons without allowing second rescatterings
        if (allowMultipleRescatterings)
          queueDecaysAndRescatterings(event, oldSize, candidates);
      }
    }

    // Third part: include Bose-Einstein effects among current particles.
    if (doBoseEinsteinNow) {
      if (!boseEinstein.shiftEvent(event)) return false;
      doBoseEinsteinNow = false;
    }
 
    // Fourth part: sequential decays also of long-lived particles.
    if (doDecay) {
      decaysProducedPartons = decays.decayAll(event);
    }

  // @TODO If not done, the next time around, rescatterings are not in order
  // Normally done first time around, but sometimes not (e.g. Upsilon).
  } while (decaysProducedPartons);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Calculate possible decays and rescatterings and add all to the queue.

void HadronLevel::queueDecaysAndRescatterings(Event& event, int iStart, 
  priority_queue<HadronLevel::PriorityNode>& queue)
{
  // @TODO Stress test and optimise if needed
  for (int iFirst = iStart; iFirst < event.size(); ++iFirst) 
  {
    Particle& hadA = event[iFirst];

    // @TODO Have the right upper bound here
    if (doDecay && hadA.isFinal() && hadA.canDecay() && hadA.mayDecay()
    && (hadA.mWidth() > widthSepBE || hadA.id() == 311))
      queue.push(PriorityNode(iFirst, hadA.vDec()));

    if (!hadA.isFinal() || !hadA.isHadron())
      continue;

    for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
      Particle& hadB = event[iSecond];
      if (!hadB.isFinal() || !hadB.isHadron())
        continue;
      
      // Continue if the two particles come from the same decay
      // @TODO: Test if this actually does something
      /*if (hadA.mother1() == hadB.mother1() 
          && (hadA.status() >= 90 && hadA.status() <= 99))
      {
        continue;
      }
*/
      
      Vec4 origin;
      if (rescatterings.calculateRescatterOrigin(iFirst, iSecond, event, origin))
        queue.push(PriorityNode(iFirst, iSecond, origin));
    }
  }
}

//--------------------------------------------------------------------------

// Allow more decays if on/off switches changed.
// Note: does not do sequential hadronization, e.g. for Upsilon.

bool HadronLevel::moreDecays( Event& event) {

  // Colour-octet onia states must be decayed to singlet + gluon.
  if (!decayOctetOnia(event)) return false;

  // Loop through all entries to find those that should decay.
  int iDec = 0;
  do {
    if ( event[iDec].isFinal() && event[iDec].canDecay()
      && event[iDec].mayDecay() ) decays.decay( iDec, event);
  } while (++iDec < event.size());

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Decay colour-octet onium states.

bool HadronLevel::decayOctetOnia(Event& event) {

  // Loop over particles and decay any onia encountered.
  for (int iDec = 0; iDec < event.size(); ++iDec)
  if (event[iDec].isFinal()
    && particleDataPtr->isOctetHadron(event[iDec].id())) {
    if (!decays.decay( iDec, event)) return false;

    // Set colour flow by hand: gluon inherits octet-onium state.
    int iGlu = event.size() - 1;
    event[iGlu].cols( event[iDec].col(), event[iDec].acol() );
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Trace colour flow in the event to form colour singlet subsystems.
// Option will keep junctions in the remainsJunction list,
// and not eliminate any junctions by insertion.

bool HadronLevel::findSinglets(Event& event, bool keepJunctions) {

  // Clear up storage.
  colConfig.clear();

  // Find a list of final partons and of all colour ends and gluons.
  if (colTrace.setupColList(event)) return true;

  // Begin arrange the partons into separate colour singlets.

  // Junctions: loop over them, and identify kind.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
  if (event.remainsJunction(iJun)) {
    if (!keepJunctions) event.remainsJunction(iJun, false);
    int kindJun = event.kindJunction(iJun);
    iParton.resize(0);

    // Loop over junction legs.
    for (int iCol = 0; iCol < 3; ++iCol) {
      int indxCol = event.colJunction(iJun, iCol);
      iParton.push_back( -(10 + 10 * iJun + iCol) );
      // Junctions: find color ends.
      if (kindJun % 2 == 1 && !colTrace.traceFromAcol(indxCol, event, iJun,
        iCol, iParton)) return false;
      // Antijunctions: find anticolor ends.
      if (kindJun % 2 == 0 && !colTrace.traceFromCol(indxCol, event, iJun,
        iCol, iParton)) return false;
    }

    // A junction may be eliminated by insert if two quarks are nearby.
    if (!keepJunctions) {
      int nJunOld = event.sizeJunction();
      if (!colConfig.insert(iParton, event)) return false;
      if (event.sizeJunction() < nJunOld) --iJun;
    }
  }

  // Open strings: pick up each colour end and trace to its anticolor end.
  while (!colTrace.colFinished()) {
    iParton.resize(0);
    if (!colTrace.traceFromCol( -1, event, -1, -1, iParton)) return false;

    // Store found open string system. Analyze its properties.
    if (!colConfig.insert(iParton, event)) return false;
  }

  // Closed strings : begin at any gluon and trace until back at it.
  while (!colTrace.finished()) {
    iParton.resize(0);
    if (!colTrace.traceInLoop(event, iParton)) return false;

    // Store found closed string system. Analyze its properties.
    if (!colConfig.insert(iParton, event)) return false;
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Extract rapidity pairs of string pieces.

vector< vector< pair<double,double> > > HadronLevel::rapidityPairs(
  Event& event) {

  // Loop over all string systems in the event.
  vector< vector< pair<double,double> > > rapPairs;
  for (int iSub = 0; iSub < int(colConfig.size()); iSub++) {
    vector< pair<double,double> > rapsNow;
    vector<int> iPartons = colConfig[iSub].iParton;

    // Special treatment for junction systems.
    if (colConfig[iSub].hasJunction) {
      // Pick smallest and largest rapidity parton.
      double ymi = 1e10;
      double yma = -1e10;
      for (int iP = 0; iP < int(iPartons.size()); iP++) {
        int iQ = iPartons[iP];
        if (iQ < 0) continue;
        if (event[iQ].id() == 21) continue;
        double yNow = yMax(event[iQ], MTINY);
        if (yNow > yma) yma = yNow;
        if (yNow < ymi) ymi = yNow;
      }
      rapsNow.push_back( make_pair(ymi, yma) );

    // Normal strings. For closed gluon loop include first-last pair.
    } else {
      int size = int(iPartons.size());
      int end  = size - (colConfig[iSub].isClosed ? 0 : 1);
      for (int iP = 0; iP < end; iP++) {
        int    i1  = iPartons[iP];
        int    i2  = iPartons[(iP+1)%size];
        double y1  = yMax(event[i1], MTINY);
        double y2  = yMax(event[i2], MTINY);
        double ymi = min(y1, y2);
        double yma = max(y1, y2);
        rapsNow.push_back( make_pair(ymi, yma) );
      }
    }
    rapPairs.push_back(rapsNow);
  }
  // Done.
  return rapPairs;
}

//==========================================================================

} // end namespace Pythia8
