#include <iostream>
#include "Pythia8/Pythia.h"


using namespace Pythia8;

int main(int argc, const char *argv[]) {
	Pythia pythia(string("../share/Pythia8/xmldoc"), false);

	pythia.readFile("mymain.cmnd");

	pythia.init();

	pythia.next();

	pythia.stat();

	return 0;
}