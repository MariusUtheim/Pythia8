// test112.cc
// Check that forceTimeShower works.
// For Ferber and Rostomyan at Belle II.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {

  // Number of events and number to list. Mass.
  int nEvent = 100;
  int nList  = 5;
  double eCM = 12.;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Switch off automatic event listing in favour of manual.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Initialize.
  pythia.init();

  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Reset event record to allow for new event.
    event.reset();

    // Store intermediate and outgoing particles.
    event.append( 23, -22, 0, 0, 2, 3,   0,   0, 0., 0.,  0., eCM, eCM);
    event.append(  2,  23, 1, 0, 0, 0, 101,   0, 0., 0.,  0.5 * eCM,
      0.5 * eCM,  0., eCM);
    event.append( -2,  23, 1, 0, 0, 0,   0, 101, 0., 0., -0.5 * eCM,
      0.5 * eCM,  0., eCM);

    // Force parton shower.
    pythia.forceTimeShower( 2, 3, eCM);

    // Generate hadronization. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // List first few events.
    if (iEvent < nList) event.list();

  // Done.
  }
  return 0;

}
