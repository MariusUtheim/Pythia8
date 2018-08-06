#include <iostream>
#include "Pythia8/Pythia.h"

using namespace Pythia8;

class MyRescatteringUserHook : public UserHooks
{

	bool canVetoRescatteringInteraction() const override { return true; }

	bool doVetoRescatteringInteraction(Particle& p1, Particle& p2, Vec4& origin) override {
		cout << "Got rescattering interaction" << endl;
		return false;
	}

};

int main(int argc, const char *argv[]) {

	MyRescatteringUserHook userHook;

	Pythia pythia(string("../share/Pythia8/xmldoc"), false);
	pythia.addUserHooksPtr(&userHook);

	pythia.readFile("mymain.cmnd");

	pythia.init();

	pythia.next();

	pythia.stat();

	return 0;
}