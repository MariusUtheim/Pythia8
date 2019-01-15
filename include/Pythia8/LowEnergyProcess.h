#ifndef Low_Energy_Process_H
#define Low_Energy_Process_H

#include "Pythia8/Event.h"

namespace Pythia8 {

class LowEnergyProcess {
public:

  virtual double sigmaTotal(int i1, int i2, const Event& event) const = 0;

  virtual bool collide(int i1, int i2, Event& event) = 0;

};

}

#endif