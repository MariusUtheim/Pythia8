
#include "../tests.h"


static vector<pair<int, int>> particlePairs = {
    { 2212, 2212 },
    { 2212, 2112 },
    { 2212, 2224 },
    { 2112, 2224 },
    { 3112, 2224 },

    { 2212, -2212 },
    { 2212, -2112 },
    { 2112, -2212 },
    { 2112, -2112 },
    { 3112, -2224 },
    { 3212, -2112 },
    { 3212, -2212 },

    { 2212, 211 },
    { 2212, 111 },
    { 2212, -211 },
    { 2112, 111 },
    { 2112, 113 },
    { 2112, 313 }
};

static vector<double> eCMs = { 0., 1., 1.5, 2., 2.5, 3., 3.5, 4.5, 6., 8., 10. };


void test_Particle_sigmas_equals_antiparticle_sigmas() {
  cout << "Running test: Particle sigmas equals antiparticle sigmas" << endl;

  Pythia pythia("../share/Pythia8/xmldoc", false);
  pythia.readString("Print:quiet = on");
  pythia.readString("NonPerturbative:all = on");
  if (!pythia.init()) {
    cout << "Pythia failed to init.\nTEST FAILED." << endl;
    return;
  }

  HadronWidths& hadronWidths = pythia.hadronWidths;
  LowEnergySigma& sigma = pythia.lowEnergySigma;

  for (auto p : particlePairs) {
    int idA = p.first, idB = p.second;
    int antiA = pythia.particleData.hasAnti(idA) ? -idA : idA;
    int antiB = pythia.particleData.hasAnti(idB) ? -idB : idB;
    for (int proc = 1; proc <= 9; ++proc)
      for (auto eCM : eCMs) {
        if (sigma.sigmaPartial(idA, idB, eCM, proc) != sigma.sigmaPartial(antiA, antiB, eCM, proc))
          printf("Particle-antiparticle symmetry broken for %d + %d @ %f GeV during process %d\n", idA, idB, eCM, proc);
      }
  }
  cout << "Test complete." << endl;
}


void test_Particles_pick_same_process_as_antiparticles() {
  cout << "Running test: Particles pick same process as antiparticles" << endl;

  Pythia pythia("../share/Pythia8/xmldoc", false);
  pythia.readString("Print:quiet = on");
  pythia.readString("NonPerturbative:all = on");
  if (!pythia.init()) {
    cout << "Pythia failed to init.\nTEST FAILED." << endl;
    return;
  }

  HadronWidths& hadronWidths = pythia.hadronWidths;
  LowEnergySigma& sigma = pythia.lowEnergySigma;

  // Sample processes for particles
  pythia.rndm.init(1);
  vector<int> particleProcs;
  for (auto p : particlePairs) {
    int idA = p.first;
    int idB = p.second;
    for (auto eCM : eCMs)
      for (int nRepeat = 0; nRepeat  < 10; nRepeat++)
        particleProcs.push_back(sigma.pickProcess(idA, idB, eCM));
  }

  // Sample processes for antiparticles
  pythia.rndm.init(1);
  vector<int> antiParticleProcs;
  for (auto p : particlePairs) {
    int idA = pythia.particleData.antiId(p.first);
    int idB = pythia.particleData.antiId(p.second);
    for (auto eCM : eCMs)
      for (int nRepeat = 0; nRepeat  < 10; nRepeat++)
        particleProcs.push_back(sigma.pickProcess(idA, idB, eCM));
  }

  // Check that particles and antiparticles picked the same processes
  int i = 0;
  for (auto p : particlePairs)
  for (auto eCM : eCMs)
  for (int nRepeat = 0; nRepeat < 10; nRepeat++) {
    if (particleProcs[i] != antiParticleProcs[i]) {
      printf("Particle-antiparticle symmetry broken for %d + %d @ %f GeV.\n"
             "Particle picked %d, antiparticle picked %d\n", 
             p.first, p.second, eCM, particleProcs[i], antiParticleProcs[i]);
    }
    ++i;
  }

  cout << "Test complete." << endl << endl;
}

