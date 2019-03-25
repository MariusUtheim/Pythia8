#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main() {

  int iDiv = 4;
  vector<int> debug;
  for (int i = 0; i <= iDiv; ++i) debug.push_back(i);
  int& debugRef = debug[iDiv];
  int* debugPtr = &debug[iDiv];
  for (int i = iDiv + 1; i < 20; ++i) {
    debug.push_back(i);
    ++debug[iDiv];
    cout << " value = " << debug[iDiv] << ", ref = " << debugRef 
         << ", ptr = " << *debugPtr << endl;
    if (i == 12) {
      vector<int> debug2;
      for (int j = 0; j <= iDiv; ++j) debug2.push_back(j + 100);
    }
  } 

  // Done.
  return 0;
}
