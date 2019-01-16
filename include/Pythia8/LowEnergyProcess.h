#ifndef Low_Energy_Process_H
#define Low_Energy_Process_H

#include "Pythia8/Event.h"

namespace Pythia8 {

class LowEnergyProcess {
public:

  virtual bool collide(int i1, int i2, Event& event) const = 0;

};

}

#endif