#include <chrono>
#include <iostream>
#include "tests.h"

void singleEvent()
{
  Pythia pythia;
  pythia.readFile("mymain.cmnd");
  pythia.init();

  pythia.next();

}

void fingerprint()
{
  Pythia pythia("../share/Pythia8/xmldoc", false);
  pythia.readFile("mymain.cmnd");
  pythia.readString("Print:quiet = on");
  pythia.init();

  pythia.next();
  cout << pythia.event[pythia.event.size() - 1].p();

/*
  while (true)
  {
    if (!pythia.next()) continue;

    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 553)  {
        cout << pythia.event[pythia.event.size() - 1].p();
      return;}
  }
*/
}

int main(int argc, const char *argv[]) {

  fingerprint();

  return 0;
}