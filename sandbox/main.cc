#include <chrono>
#include <iostream>
#include "Pythia8/Pythia.h"

using namespace Pythia8;
using namespace std::chrono;



void main_directionality() {

	Pythia pythia(string("../share/Pythia8/xmldoc"), false);

	pythia.readFile("mymain.cmnd");
	pythia.init();

	int nOutwards = 0, nInwards = 0;
	double eps = 1.0e-12;
	Hist cosalpha("cosalpha", 100, -1.0 - eps, 1.0 + eps);
	double rmax = pythia.settings.parm("Rescattering:maxRadius");

	for (int iEvent = 0; iEvent < 1000; ++iEvent)
	{
		pythia.next();

		for (int i = 0; i < pythia.event.size(); ++i)
		{
			Particle p = pythia.event[i];

			if (p.isFinal() && p.isHadron() && p.vProd().pAbs() < rmax
					&& p.vProd().pAbs() > 0)
			{
				double c = dot3(p.vProd() / p.vProd().pAbs(), p.p() / p.p().pAbs());
				
				if (c >= 0)
					++nOutwards;
				else
					++nInwards;

				cosalpha.fill(c);
			}
		}
	}

	cout << cosalpha;
}


void main_timeRescattering() {

	int nRuns = 10000;


	auto tBefore = high_resolution_clock::now();
	auto tAfter = high_resolution_clock::now();
	duration_cast<milliseconds>(tAfter - tBefore).count();

	int basicTime, rescTime, rerescTime;


	Pythia pythiaB(string("../share/Pythia8/xmldoc"), false);
	pythiaB.readFile("mymain.cmnd");
	pythiaB.init();

	tBefore = high_resolution_clock::now();
	for (int i = 0; i < nRuns; ++i) 
		pythiaB.next();
	tAfter = high_resolution_clock::now();
	basicTime = duration_cast<milliseconds>(tAfter - tBefore).count();


	Pythia pythiaR(string("../share/Pythia8/xmldoc"), false);
	pythiaR.readFile("mymain.cmnd");
	pythiaR.readString("Rescattering:rescattering = on");
	pythiaR.readString("Rescattering:allowSecondRescattering = off");
	pythiaR.init();

	tBefore = high_resolution_clock::now();
	for (int i = 0; i < nRuns; ++i) 
		pythiaR.next();
	tAfter = high_resolution_clock::now();
	rescTime = duration_cast<milliseconds>(tAfter - tBefore).count();


	Pythia pythia2(string("../share/Pythia8/xmldoc"), false);
	pythia2.readFile("mymain.cmnd");
	pythia2.readString("Rescattering:rescattering = on");
	pythia2.readString("Rescattering:allowSecondRescattering = on");
	pythia2.init();

	tBefore = high_resolution_clock::now();
	for (int i = 0; i < nRuns; ++i) 
		pythia2.next();
	tAfter = high_resolution_clock::now();
	rerescTime = duration_cast<milliseconds>(tAfter - tBefore).count();


	cout << endl << " ==== RESULTS ==== " << endl
			 << " Basic " << basicTime << endl
			 << " Resca " << rescTime << endl 
			 << " ReRes " << rerescTime << endl
			 << endl ;


}

void main_doTest() {
	Pythia pythia(string("../share/Pythia8/xmldoc"), false);

	pythia.readFile("mymain.cmnd");
	pythia.init();

	for (int i = 0; i < 100; ++i)
		pythia.next();

}


int main(int argc, const char *argv[]) {

	main_doTest();

	return 0;
}