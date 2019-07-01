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

using std::priority_queue;

//--------------------------------------------------------------------------

// Node for ordering scatterings and decays

class HadronLevel::PriorityNode {
public:
  PriorityNode(int iDecayIn, Vec4 originIn)
    : i1(iDecayIn), i2(0), origin(originIn) {}
  PriorityNode(int i1In, int i2In, Vec4 originIn)
    : i1(i1In), i2(i2In), origin(originIn) {}
  
  bool isDecay() { return i2 == 0; }

  // Priority comparison to be used by priority_queue
  // Note that lower t means higher priority!
  // @FUTURE: allow the user to pick the comparer (time-ordering problem)?
  bool operator<(const PriorityNode& r) const
  { return origin.e() > r.origin.e(); }

  int i1, i2;
  Vec4 origin;
};

//--------------------------------------------------------------------------

// Find settings. Initialize HadronLevel classes as required.

bool HadronLevel::init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
  ParticleData* particleDataPtrIn, HadronWidths* hadronWidthsPtrIn,
  Couplings* couplingsPtrIn, TimeShower* timesDecPtr,
  RHadrons* rHadronsPtrIn, DecayHandler* decayHandlePtr,
  vector<int> handledParticles, LowEnergySigma* lowEnergySigmaPtrIn,
  UserHooks* userHooksPtrIn) {

  // Save pointers.
  infoPtr           = infoPtrIn;
  particleDataPtr   = particleDataPtrIn;
  rndmPtr           = rndmPtrIn;
  couplingsPtr      = couplingsPtrIn;
  rHadronsPtr       = rHadronsPtrIn;
  userHooksPtr      = userHooksPtrIn;
  lowEnergySigmaPtr = lowEnergySigmaPtrIn;

  // Main flags.
  doHadronize     = settings.flag("HadronLevel:Hadronize");
  doDecay         = settings.flag("HadronLevel:Decay");
  doRescatter     = settings.flag("HadronLevel:Rescatter");
  doBoseEinstein  = settings.flag("HadronLevel:BoseEinstein");

  // Check settings for rescattering
  if (doRescatter && !settings.flag("Fragmentation:setVertices")) {
    infoPtr->errorMsg("Error in HadronLevel::init: HadronLevel:Rescatter "
      "is on, but Fragmentation:setVertices is off");
    return false;
  }

  // Boundary mass between string and ministring handling.
  mStringMin      = settings.parm("HadronLevel:mStringMin");

  // For junction processing.
  eNormJunction   = settings.parm("StringFragmentation:eNormJunction");

  // Allow R-hadron formation.
  allowRH         = settings.flag("RHadrons:allow");

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

  // Allow decayed/rescattered particles to rescatter again
  scatterManyTimes = settings.flag("Rescattering:scatterManyTimes");

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

  // Initialize particle decays.
  decays.init(infoPtr, settings, particleDataPtr, rndmPtr,
    hadronWidthsPtrIn, couplingsPtr, timesDecPtr, &flavSel,
    decayHandlePtr, handledParticles);

  // Initialize low energy.
  lowEnergyProcess.init(infoPtr, settings, rndmPtr, 
    particleDataPtr, hadronWidthsPtrIn, &stringFrag, &ministringFrag);

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

// Calculate possible decays and rescatterings and add all to the queue.

void HadronLevel::queueDecResc(Event& event, int iStart, 
  priority_queue<HadronLevel::PriorityNode>& queue)
{
  for (int iFirst = iStart; iFirst < event.size(); ++iFirst) 
  {
    Particle& hadA = event[iFirst];
    if (!hadA.isFinal() || !hadA.isHadron())
      continue;
    
    // Queue hadrons that should decay
    if (doDecay && hadA.canDecay() && hadA.mayDecay()
    && (hadA.mWidth() > widthSepBE || hadA.id() == 311)) 
      queue.push(PriorityNode(iFirst, hadA.vDec()));

    // Loop over particle pairs
    for (int iSecond = 0; iSecond < iFirst; ++iSecond) {
      
      // @4TS verify the logic here
      Particle& hadB = event[iSecond];
      if (!hadB.isFinal() || !hadB.isHadron())
        continue;

      // Abort early if particles are moving away from each other
      if (dot3(hadB.p() / hadB.e() - hadA.p() / hadA.e(), 
               hadB.vProd() - hadA.vProd() ) > 0)
        continue;

      // Set up positions for each particle in their CM frame
      RotBstMatrix frame;
      frame.toCMframe(hadA.p(), hadB.p());

      Vec4 vA = hadA.vProd(), vB = hadB.vProd();
      Vec4 pA = hadA.p(),     pB = hadB.p();

      vA.rotbst(frame); vB.rotbst(frame);
      pA.rotbst(frame); pB.rotbst(frame);

      // Offset particles to position when the last particle is created
      double t0 = max(vA.e(), vB.e());
      double zA = vA.pz() + (t0 - vA.e()) * pA.pz() / pA.e();
      double zB = vB.pz() + (t0 - vB.e()) * pB.pz() / pB.e();

      // Abort if particles have already passed each other
      if (zA >= zB)
        continue;

      // Calculate sigma and abort if impact parameter is too large
      double eCM = (pA + pB).mCalc();
      double sigma = lowEnergySigmaPtr->sigmaTotal(hadA.id(), hadB.id(), eCM);
      if ((vA - vB).pT2() > MB2MMSQ * sigma / M_PI)
        continue;

      // Calculate collision origin and transform to lab frame
      double tColl = t0 - (zB - zA) / (pB.pz() / pB.e() - pA.pz() / pA.e());
      Vec4 origin(0.5 * (vA.px() + vB.px()), 0.5 * (vA.py() + vB.py()),
                  zA + pA.pz() / pA.e() * (tColl - t0), tColl);
      frame.invert();
      origin.rotbst(frame);

      // Queue hadron that should rescatter
      queue.push(PriorityNode(iFirst, iSecond, origin));
    }
  }
}

//--------------------------------------------------------------------------

// Hadronize and decay the next parton-level.

bool HadronLevel::next(Event& event) {

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
  bool doBoseEinsteinNow = doBoseEinstein;
  bool decaysCausedHadronization;
  do {
    decaysCausedHadronization = false;

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
      decaysCausedHadronization = decays.decayAll(event, widthSepBE);
    }
    // If rescattering is on, decays/rescatterings must happen in order
    else if (doRescatter) {

      priority_queue<PriorityNode> candidates; 

      queueDecResc(event, 0, candidates);

      while (!candidates.empty()) {
        PriorityNode node = candidates.top();
        candidates.pop();

        // Abort if either particle has already interacted elsewhere
        if (!event[node.i1].isFinal() 
        || (!node.isDecay() && !event[node.i2].isFinal()))
          continue;

        int oldSize = event.size();

        // Perform the queued action
        if (node.isDecay()) {
          decays.decay(node.i1, event);
          // @TBD If there is moreToDo, those things should also be handled in order?
          if (decays.moreToDo()) decaysCausedHadronization = true;
        }
        else {
          double eCM = (event[node.i1].p() + event[node.i2].p()).mCalc();
          int process = lowEnergySigmaPtr->pickProcess(event[node.i1].id(),
            event[node.i2].id(), eCM);
          if (process != 0)
            lowEnergyProcess.collide(node.i1, node.i2,
                                    process, event, node.origin);
        }

        // Check for new interactions
        if (scatterManyTimes)
          queueDecResc(event, oldSize, candidates);
        // @4TS Verify logic, @TODO test that the right hadrons decay before BE
        else if (doDecay) {
          // If multiple rescattering is off, particles can still decay again
          for (int i = oldSize; i < event.size(); ++i) {
            if (event[i].isFinal() && event[i].isHadron()
            && event[i].canDecay() && event[i].mayDecay()
            && (event[i].mWidth() > widthSepBE || event[i].id() == 311)) {
              decays.decay(i, event);
              if (decays.moreToDo())
                decaysCausedHadronization = true;
            }
          }
        }
      }
    }

    // Third part: include Bose-Einstein effects among current particles.
    if (doBoseEinsteinNow) {
      if (!boseEinstein.shiftEvent(event)) return false;
      doBoseEinsteinNow = false;
    }
 
    // Fourth part: sequential decays also of long-lived particles.
    if (doDecay) {
      if (decays.decayAll(event))
        decaysCausedHadronization = true;
    }

  // @TBD Do we need special considerations for dealing with more rescatters?
  // Normally done first time around, but sometimes not.
  // (e.g. Upsilon decay can cause create unstable hadrons).
  } while (decaysCausedHadronization);

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Allow more decays if on/off switches changed.
// Note: does not do sequential hadronization, e.g. for Upsilon.

bool HadronLevel::moreDecays( Event& event) {

  // @TDB: It so happens that HadronLevel::moreDecays is called from 
  //   Pythia::nextNonPert, and there is no other check to whether doDecays 
  //   is set. Perhaps there should be such a case somewhere, or perhaps 
  //   this function should not actually be called at that point.
  if (!doDecay)
    return true;

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
