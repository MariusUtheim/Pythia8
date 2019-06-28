#ifndef MYTESTS_H
#define MYTESTS_H

#include "Pythia8/Pythia.h"

using namespace Pythia8;

const vector<int> typeCodes =    { 211,   321,  313,   2212, 333,   3312,   3224,      3324    };
const vector<string> typeNames = { "pi+", "K+", "K*0", "p",  "phi", "Chi-", "Sigma*+", "Chi*0" };
const vector<double> datapT =    { 0.47,  0.77, 1.01,  0.90, 1.07,  1.21,   1.15,      1.31    };  

void test_DeltapT_per_collision();
void test_compareDeltapT();
void test_pT_distributions();
void test_ClebschGordan();

void test_Particle_sigmas_equals_antiparticle_sigmas();
void test_Particles_pick_same_process_as_antiparticles();

#endif
