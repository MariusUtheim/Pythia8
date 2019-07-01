
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergySigma.h"

namespace Pythia8 {

//==========================================================================

// The LowEnergySigma class.

//--------------------------------------------------------------------------

void LowEnergySigma::init(Info* infoPtrIn, Settings& settings, Rndm* rndmPtrIn,
    ParticleData* particleDataPtrIn, HadronWidths* hadronWidthsPtrIn) {
    
  // Store pointers
  infoPtr           = infoPtrIn;
  particleDataPtr   = particleDataPtrIn; 
  rndmPtr           = rndmPtrIn;
  hadronWidthsPtr = hadronWidthsPtrIn;

  // @TODO: Maybe just copy the relevant parts from sigmaSaSDL
  sigmaSaSDL.init(infoPtrIn, settings, particleDataPtrIn, nullptr);

  // Initialize map of resonance particles
  for (int id : hadronWidthsPtr->getResonances()) {

    // Insert id in signature map
    int signature = getSignature(particleDataPtr->isHadron(id), 
      particleDataPtr->chargeType(id), particleDataPtr->nQuarksInCode(id, 3));
    auto iter = signatureToParticles.find(signature);
    if (iter != signatureToParticles.end())
      iter->second.push_back(id);
    else
      signatureToParticles.emplace(signature, vector<int>{id});
  }

}

//--------------------------------------------------------------------------

// Returns int representing the overall process type:
//  0: Collision not implemented (ids will not be ordered)
//  1: BB
//  2: BBbar
//  3: XM
//  4: XbarM
// 
// The canonical ordering of A and B satisfies two criteria:
//   1) |A| >= |B|, and 2) A > 0
// 
// Implications:
//  - In BB, the antiparticle cross sections are the same as the particle ones,
//    so BbarBbar is replaced by BB
//  - In BBbar, A is always the particle and B is the antiparticle. Again, 
//    signs are flipped if necessary
//  - In XM, X is a baryon or meson and B is always a meson.
//
int LowEnergySigma::canonicalForm(int& idA, int& idB) const {
  if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB))
    return 0;
  // Explicitly don't handle c or b hadrons
  if (abs(particleDataPtr->heaviestQuark(idA)) > 3 
   || abs(particleDataPtr->heaviestQuark(idB)) > 3)
    return 0;

  // Ensure |A| >= |B|
  if (abs(idA) < abs(idB))
    swap(idA, idB);

  // Ensure A > 0
  bool flipSign = idA < 0;
  if (flipSign) { 
    idA = -idA;
    idB = particleDataPtr->antiId(idB);
  }

  // Get id of overall collision type
  if (particleDataPtr->isMeson(idB))
    return flipSign ? 4 : 3; // XM
  else if (idB < 0)
    return 2; // BBbar
  else
    return 1; // BB
}

//--------------------------------------------------------------------------

double LowEnergySigma::sigmaTotal(int idA, int idB, double eCM) const {
 
  // If energy is less than the hadron masses, return 0.
  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return 0.;

  switch (canonicalForm(idA, idB)) {
    case 0: return 0.; // E.g. if particles are not hadrons
    case 1: return BBTotal(idA, idB, eCM);
    case 2: return BBbarTotal(idA, idB, eCM);
    case 3: case 4: return XMTotal(idA, idB, eCM);
    
    // This should never occur
    default: 
      infoPtr->errorMsg("Warning in LowEnergySigma::sigmaTotal: "
        "Unhandled case encountered. This indicates a bug in LowEnergySigma.");
      return 0.;
  }
}

//--------------------------------------------------------------------------

bool LowEnergySigma::sigmaPartial(int idA, int idB, double eCM,
  vector<int>& procsOut, vector<double>& sigmasOut) const {

  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return false;

  int process = canonicalForm(idA, idB);
  switch (process) {
    case 0: // E.g. if particles are not hadrons
      return false;

    case 1: // BB
      procsOut = { 1, 2, 3, 4, 5, 7 };
      sigmasOut = { BBNonDiff(idA, idB, eCM), BBElastic(idA, idB, eCM),
        BBDiffractiveAX(idA, idB, eCM), BBDiffractiveXB(idA, idB, eCM),
        BBDiffractiveXX(idA, idB, eCM), BBExcite(idA, idB, eCM)
      };
      return true;

    case 2: // BBbar
      procsOut = { 1, 2, 3, 4, 5, 8 };
      sigmasOut = { BBbarNonDiff(idA, idB, eCM), BBbarElastic(idA, idB, eCM),
        BBbarDiffractiveAX(idA, idB, eCM), BBbarDiffractiveXB(idA, idB, eCM),
        BBbarDiffractiveXX(idA, idB, eCM), BBbarAnnihilation(idA, idB, eCM)
      };
      return true;

    case 3: case 4: // XM
      procsOut = { 1, 2 };
      sigmasOut = { XMNonDiffractive(idA, idB, eCM), XMElastic(idA, idB, eCM)};
      // Add resonances with nonzero cross sections
      for (auto idR : possibleResonances(idA, idB)) {
        double sigmaRes = XMResonantPartial(idA, idB, idR, eCM);
        if (sigmaRes != 0.) {
          if (process == 4) idR = particleDataPtr->antiId(idR);
          procsOut.push_back(idR);
          sigmasOut.push_back(sigmaRes);
        }
      }
      return true;

    default: 
      infoPtr->errorMsg("Warning in LowEnergySigma::sigmaPartial: "
        "Unhandled case encountered. This indicates a bug in LowEnergySigma.");
      return false;
  }
}

//--------------------------------------------------------------------------

double LowEnergySigma::sigmaPartial(int idA, int idB, double eCM,
  int type) const {

  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return 0.;

  switch (canonicalForm(idA, idB)) {
    case 0: // E.g. if particles are not hadrons
      return 0.;

    case 1: // BB
      switch (type) {
        case 0: return BBTotal(idA, idB, eCM);
        case 1: return BBNonDiff(idA, idB, eCM);
        case 2: return BBElastic(idA, idB, eCM);
        case 3: return BBDiffractiveAX(idA, idB, eCM);
        case 4: return BBDiffractiveXB(idA, idB, eCM);
        case 5: return BBDiffractiveXX(idA, idB, eCM);
        case 7: return BBExcite(idA, idB, eCM);
        default: return 0.;
      }

    case 2: // BBbar
      switch (type) {
        case 0: return BBbarTotal(idA, idB, eCM);
        case 1: return BBbarNonDiff(idA, idB, eCM);
        case 2: return BBbarElastic(idA, idB, eCM);
        case 3: return BBbarDiffractiveAX(idA, idB, eCM);
        case 4: return BBbarDiffractiveXB(idA, idB, eCM);
        case 5: return BBbarDiffractiveXX(idA, idB, eCM);
        case 8: return BBbarAnnihilation(idA, idB, eCM);
        default: return 0.;
      }
    
    case 3: case 4: // XM
      switch (type) {
        case 1: return XMNonDiffractive(idA, idB, eCM);
        case 2: return XMElastic(idA, idB, eCM);
        case 9: return XMResonant(idA, idB, eCM);
        default:
          return abs(type) > 100 ? XMResonantPartial(idA, idB, type, eCM) : 0.;
      }

    default: 
      infoPtr->errorMsg("Warning in LowEnergySigma::sigmaPartial: "
        "Unhandled case encountered. This indicates a bug in LowEnergySigma.");
      return 0.;
  }
}

//--------------------------------------------------------------------------

int LowEnergySigma::pickProcess(int idA, int idB, double eCM) {
  vector<int> processes;
  vector<double> sigmas;
  if (!sigmaPartial(idA, idB, eCM, processes, sigmas))
    return 0;

  return processes[rndmPtr->pick(sigmas)];
}

//--------------------------------------------------------------------------

int LowEnergySigma::pickResonance(int idA, int idB, double eCM) {
  int process = canonicalForm(idA, idB);
  if (process == 0 || process == 1 || process == 2) 
    return 0;
  
  // Get resonances with relative frequencies
  vector<int> resonances = possibleResonances(idA, idB);
  vector<double> sigmas(resonances.size());
  
  bool foundAny = false;
  for (size_t i = 0; i < resonances.size(); ++i) {
    double sigmaRes = XMResonantPartial(idA, idB, resonances[i], eCM);
    if (sigmaRes != 0) foundAny = true;
    sigmas[i] = sigmaRes;
  }
  
  // If no resonances are available, return 0
  if (!foundAny)
    return 0;

  // Pick resonance at random
  int resPick = resonances[rndmPtr->pick(sigmas)];
  // Change to antiparticle if the canonical ordering changed signs
  return (process == 3) ? resPick : particleDataPtr->antiId(resPick);
}

//--------------------------------------------------------------------------

// Helper functions

double LowEnergySigma::aqm(int idA, int idB) const {
  double mesA = particleDataPtr->isMeson(idA);
  double mesB = particleDataPtr->isMeson(idB);
  return 40 * pow(2./3., mesA + mesB)
    * (1 - 0.4 * abs(particleDataPtr->nQuarksInCode(idA, 3)) / (mesA ? 2 : 3))
    * (1 - 0.4 * abs(particleDataPtr->nQuarksInCode(idB, 3)) / (mesB ? 2 : 3));
}

double LowEnergySigma::aqmNN() const {
  return 40;
}

static double ReggeFit(double z, double y1, double y2, double s) {
  return z + 0.308 * pow2(log(s / 28.998)) 
       + y1 * pow(s, -0.458) - y2 * pow(s, -0.545);
}

static double HERAFit(double a, double b, double n, double c, double d, double p) {
  return a + b * pow(p, n) + c * pow2(log(p)) + d * log(p);
}

//--------------------------------------------------------------------------

// Baryon-Baryon section

/** Todo list for BB:
 * @TODO Check that the tables are correct and sufficiently smooth, compare w/ data
 * @TODO Implement strangeness exchange
 * @TODO Implement parametrisation for Lambda+p and Sigma+p special cases
 * @TODO Implement D+N and D+D collisions
 * @TODO Do something abour charm and bottom? 
 **/

// === Begin interpolation data for NN cross sections ===

static Interpolator ppTotalData(1.88, 5.0, {
  314.914, 60.7018, 30.4889, 25.1787, 24.5172, 24.125, 22.744, 24.151, 
  24.4887, 25.5162, 26.0655, 28.7429, 31.5493, 35.0329, 38.0288, 
  42.1204, 43.5899, 45.5891, 46.9695, 47.4814, 47.0977, 48.031, 
  47.7431, 46.9851, 47.6045, 47.4623, 47.0786, 47.2336, 47.8374, 
  46.3766, 46.6738, 46.9709, 46.8992, 46.5721, 46.1646, 45.6672, 
  45.4572, 45.4466, 45.4359, 44.6773, 43.4844, 44.1945, 44.9942, 
  45.0851, 44.5214, 44.373, 44.2346, 44.3895, 44.95, 44.7822, 42.7161,
  42.824, 43.0213, 42.8924, 42.6889, 42.6083, 42.5568, 42.5052, 
  42.4028, 42.2919, 42.1811, 41.702, 41.1733, 41.1347, 41.4545, 
  41.7743, 41.8671, 41.7377, 41.6083, 41.4789, 41.3854, 41.3247, 
  41.2641, 41.2035, 41.1429, 41.0823, 41.0217, 40.961, 40.9004, 
  40.8398, 40.7792, 40.7186, 40.658, 40.5973, 40.5367, 40.4761, 
  40.4155, 40.3549, 40.2942, 40.2336, 40.173, 40.1124, 40.0518, 
  39.9912, 39.9305, 39.8699, 39.8093, 39.7487, 39.6881, 39.6275,
  39.5668, 39.5062, 39.4456, 39.385, 39.3244, 39.2638, 39.2031, 
  39.1425, 39.0819, 39.0213, 38.9607, 38.9, 38.8394, 38.7788, 38.7182, 
  38.6576, 38.597, 38.5363, 38.4757, 38.4151, 38.3545, 38.2939, 
  38.2333, 38.1726, 38.112, 38.0514, 37.9908, 37.9302, 37.8696, 
  37.8089, 37.7483, 37.6877, 37.6271, 37.5665, 37.5058, 37.4452, 
  37.3846, 37.324, 37.2634, 37.2028, 37.1421, 37.0815, 37.0209, 
  36.9603, 36.8997, 36.8391, 36.7784, 36.7178, 36.6572, 36.5966, 
  36.536, 36.4754, 36.4147, 36.3541, 36.2935, 36.2329, 36.1723
});

static Interpolator pnTotalData(1.88, 5.0, {
  1135.35, 170.185, 84.5462, 58.9083, 47.5717, 43.3222, 40.1235, 
  37.2117, 34.7469, 34.4045, 34.0621, 34.067, 34.6442, 35.2214, 
  35.8313, 36.8508, 37.8704, 38.6998, 38.9281, 39.1564, 39.3846, 
  39.7081, 40.1064, 40.475, 40.8079, 41.1407, 41.4218, 41.5763, 
  41.7324, 41.8934, 42.0543, 42.2153, 42.3605, 42.4568, 42.5531, 
  42.6495, 42.7458, 42.8188, 42.8745, 42.9302, 42.9702, 43.0034, 
  43.0382, 43.0772, 43.1136, 43.1069, 43.0772, 43.0018, 42.9265, 
  42.8511, 42.802, 42.7533, 42.7046, 42.6559, 42.6072, 42.5585, 
  42.5173, 42.4955, 42.4738, 42.452, 42.4302, 42.4084, 42.3866, 42.361, 
  42.3353, 42.3096, 42.2838, 42.2581, 42.2323, 42.2066, 42.1809, 
  42.1551, 42.1294, 42.1037, 42.0779, 42.0522, 42.0334, 42.0186, 
  42.0037, 41.9888, 41.974, 41.9591, 41.9442, 41.9293, 41.9145, 
  41.8996, 41.8847, 41.8698, 41.855, 41.8401, 41.8252, 41.8103, 
  41.7955, 41.7806, 41.7657, 41.7508, 41.736, 41.7211, 41.7062, 
  41.6913, 41.6765, 41.6616, 41.6467, 41.6318, 41.617, 41.6021, 
  41.5872, 41.5724, 41.5575, 41.5426, 41.5277, 41.5129, 41.498, 
  41.4831, 41.4682, 41.4534, 41.4385, 41.4236, 41.4087, 41.3939, 
  41.379, 41.3641, 41.3492, 41.3344, 41.3195, 41.3046, 41.2897, 
  41.2749, 41.26, 41.2451, 41.2302, 41.2154, 41.2005, 41.1856, 41.1708,
  41.1559, 41.141, 41.1261, 41.1113, 41.0964, 41.0815, 41.0666, 
  41.0518, 41.0369, 41.022, 41.0071, 40.9923, 40.9774, 40.9625, 
  40.9476, 40.9328, 40.9179, 40.903, 40.8881, 40.8733, 40.8584, 40.8435
});

static Interpolator ppElasticData(2.0, 5.0, {
  22.01, 22.7518, 23.4936, 24.2354, 24.1411, 23.9591, 24.3721, 
  24.7462, 24.8487, 24.2555, 23.803, 23.4271, 23.2217, 24.1612, 
  25.1006, 25.7232, 25.5535, 23.2231, 21.3995, 21.3974, 21.3952, 
  21.393, 21.3908, 21.3887, 21.3596, 21.0491, 20.7386, 20.4282, 
  20.1177, 19.8072, 19.4967, 19.1862, 18.8757, 18.5652, 18.2547, 
  17.9442, 17.6337, 17.4215, 17.2347, 17.0479, 16.8611, 16.6743, 
  16.4875, 16.3006, 16.1138, 15.927, 15.7297, 15.5232, 15.3168, 
  15.1103, 14.9038, 14.6973, 14.4909, 14.2844, 14.0779, 13.8715,
  13.665, 13.4585, 13.252, 13.0953, 12.9736, 12.8519, 12.7302, 12.6085,
  12.4868, 12.3652, 12.2435, 12.1218, 12.0001, 11.8784, 11.7567,
  11.635, 11.5439, 11.4965, 11.4491, 11.4017, 11.3543, 11.3069, 
  11.2595, 11.2121, 11.1648, 11.1174, 11.07, 11.0226, 11.0005, 11.0343,
  11.0681, 11.1019, 11.1357, 11.1695, 11.2033, 11.2371, 11.0865, 
  10.2933, 9.50015, 9.21152, 9.2889, 9.36629, 9.44368, 9.52106, 
  9.59845, 9.67584, 9.75323, 9.83061, 9.908, 9.98539, 10.0628, 10.1402, 
  10.2175, 10.2949, 10.3723, 10.4014, 10.3658, 10.3301, 10.2944, 
  10.2587, 10.1975, 10.1233, 10.049, 9.9747, 9.90042, 9.82615, 9.8552, 
  9.94039, 10.0256, 10.1108, 10.196, 10.2812, 10.3663, 10.4515, 
  10.4507, 10.4067, 10.3628, 10.3188, 10.2748, 10.2308, 10.1869, 
  10.1429, 10.0989, 10.0549, 10.0109, 9.96697, 9.92299, 9.87901, 
  9.95995, 10.0594, 10.1589, 10.2583, 10.3578, 10.4572, 10.5567
});

static Interpolator pnElasticData(2.0, 4.0, {
  35.9272, 27.5224, 23.087, 18.528, 14.125, 15.4334, 14.4921, 12.8764,
  12.0672, 10.7543, 8.17421
});

static Interpolator NNExciteData(3.8, 14.3, {
  28.623, 28.4801, 28.2446, 27.9334, 27.566, 27.1519, 26.709, 26.2366, 
  25.7466, 25.2436, 24.7344, 24.2218, 23.7088, 23.198, 22.6918, 
  22.1911, 21.698, 21.2129, 20.7364, 20.2698, 19.8131, 19.3665, 
  18.9308, 18.5055, 18.0907, 17.6866, 17.293, 16.9094, 16.5365, 
  16.1737, 15.8204, 15.4772, 15.1432, 14.8184, 14.5037, 14.1978, 
  13.933, 13.6465, 13.3713, 13.1129, 13.0205, 12.7895, 12.677, 12.4801, 
  12.2832, 12.0863, 11.8894, 11.6925, 11.4956, 11.2987, 11.1018, 
  10.9048, 10.7079, 10.511, 10.3141, 10.1172, 9.92031, 9.7234, 9.52649, 
  9.32958, 9.13268, 8.93577, 8.73886, 8.54195, 8.34504, 8.14814, 
  7.95123, 7.75432, 7.55741, 7.36051, 7.1636, 6.96669, 6.76978, 
  6.57287, 6.37597, 6.17906, 5.98215, 5.78524, 5.58833, 5.39143, 
  5.19452, 4.99761, 4.8007, 4.60379, 4.40689, 4.20998, 4.01307, 
  3.81616, 3.61925, 3.42235, 3.22544, 3.02853, 2.83162, 2.63471, 
  2.43781, 2.2409, 2.04399, 1.84708, 1.65018, 1.45327, 1.25636, 
  1.05945, 0.862544, 0.665636, 0.468728, 0.27182
});

// === End interpolation data for NN cross sections ===

// All diffractive processes below threshold are excitations
constexpr double NNExciteThreshold = 3.8;

// Diffractive below threshold fall linearly to zero at excitation threshold
constexpr double SaSDLThreshold = 8.;


double LowEnergySigma::BBTotal(int idA, int idB, double eCM) const {
  // Use parametrisation for pp/nn
  if ((idA == 2212 && idB == 2212) || (idA == 2112 && idB == 2112)) {
    double t = clamp((eCM - 3.) / (5. - 3.), 0., 1.);
    return (1 - t) * ppTotalData(eCM) 
               + t * ReggeFit(35.45, 42.53, 33.34, eCM * eCM);
  }
  // Use parametrisation for pn
  else if (idA == 2212 && idB == 2112)
  {
    double t = clamp((eCM - 3.) / (5. - 3.), 0., 1.);
    return (1 - t) * pnTotalData(eCM) 
               + t * ReggeFit(35.80, 40.15, 30.00, eCM * eCM);
  }
  // @TODO: Something special for Delta1232+N or Delta1232+Delta1232
  // ...
  // Use AQM + strangeness exchange for all others
  else {
    // @TODO: Add strangeness exchange? 
    return aqm(idA, idB);
  }
}

double LowEnergySigma::BBElastic(int idA, int idB, double eCM) const {
  // Fit pp/nn/pn
  if ((idA == 2112 || idA == 2212) && (idB == 2112 || idB == 2212)) {
    if (eCM < 2 * particleDataPtr->m0(2212) + particleDataPtr->m0(211))
      return BBTotal(idA, idB, eCM);

    double t = clamp((eCM - 3.) / (5. - 3.), 0., 1.);

    double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
    double s = eCM * eCM;
    double pLab = sqrt((s - pow2(mA + mB)) * (s - pow2(mA - mB))) / (2. * mB);
    
    // HERA fit is the same for others (pp and pn are simlar at high energies)
    double sigmaHERA = HERAFit(11.9, 26.9, -1.21, 0.169, -1.85, pLab);

    // Data fit at low energies is different for pp/nn or pn
    double sigmaData = (idA == idB) ? ppElasticData(eCM)  // pp/nn
                                    : pnElasticData(eCM); // pn

    return (1 - t) * sigmaData + t * sigmaHERA;
  }
  else {
    // For processes other than NN, use AQM
    return 0.039 * pow(aqm(idA, idB), 2./3.);
  }
}

double LowEnergySigma::BBNonDiff(int idA, int idB, double eCM) const {
  return BBTotal(idA, idB, eCM) - BBElastic(idA, idB, eCM)
        - BBDiffractiveAX(idA, idB, eCM) - BBDiffractiveXB(idA, idB, eCM)
        - BBDiffractiveXX(idA, idB, eCM) - BBExcite(idA, idB, eCM);
}

double LowEnergySigma::BBDiffractiveAX(int idA, int idB, double eCM) const {

  // @TODO: Diffractive scattering is only implemented for NN collisions
  if ((idA != 2212 && idA != 2112) || (idB != 2212 && idB != 2112))
    return 0.;

  // No continuous diffraction in the resonance region
  if (eCM < NNExciteThreshold) 
    return 0.;

  // Interpolate between resonance region and beginning of SaSDL region
  double t = 1.;
  if (eCM < SaSDLThreshold) {
    t = (eCM - NNExciteThreshold) / (SaSDLThreshold - NNExciteThreshold);
    eCM = SaSDLThreshold;
  }
  
  // Calculate cross section from SaSDL
  sigmaSaSDL.calcDiff(idA, idB, eCM * eCM, 
                      particleDataPtr->m0(idA), particleDataPtr->m0(idB));
  return t * sigmaSaSDL.sigAX;
}

double LowEnergySigma::BBDiffractiveXB(int idA, int idB, double eCM) const {
  return BBDiffractiveAX(idB, idA, eCM);
}

double LowEnergySigma::BBDiffractiveXX(int idA, int idB, double eCM) const {
  // @TODO: Diffractive scattering is only implemented for NN collisions
  if ((idA != 2212 && idA != 2112) || (idB != 2212 && idB != 2112))
    return 0.;

  // No continuous diffraction in the resonance region
  if (eCM < NNExciteThreshold) 
    return 0.;

  // Interpolate between resonance region and beginning of SaSDL region
  double t = 1.;
  if (eCM < SaSDLThreshold) {
    t = (eCM - NNExciteThreshold) / (SaSDLThreshold - NNExciteThreshold);
    eCM = SaSDLThreshold;
  }
  
  // Calculate cross section from SaSDL
  sigmaSaSDL.calcDiff(idA, idB, eCM * eCM, 
                      particleDataPtr->m0(idA), particleDataPtr->m0(idB));
  return t * sigmaSaSDL.sigXX;
}

double LowEnergySigma::BBExcite(int idA, int idB, double eCM) const {
  // @TODO: Excitations are only implemented for NN collisions
  if ((idA != 2212 && idA != 2112) || (idB != 2212 && idB != 2112))
    return 0.;

  // Below excitation threshold, all non-elastic processes are excitations
  if (eCM < NNExciteThreshold)
    return BBTotal(idA, idB, eCM) - BBElastic(idA, idB, eCM);

  // Above excitation threshold, parametrise excitation cross section
  if (eCM < NNExciteData.right())
    return NNExciteData(eCM);
  else
    return 0.;
}

//--------------------------------------------------------------------------

// Baryon-Antibaryon section

/**TODO list for BBbar:
 * @TODO Check that sNN is correct
 * @TODO UrQMD actually uses Regge fit instead of HERA fit for sigmaTotNN
 * @TODO sigmaTotNN and sigmaElNN do not match data well for pLab < 0.3
 * @TODO Should there be a different parametrisation for npbar?
 * @TODO Check that the aqmFactor is correct (should we use the elastic one?)
 * @TODO Figure out if the annihilation cross section parametrisation is up to date
 * @TODO Decide how to split the diffractive cross sections among specific processes
 * @TODO Compare with data cases: ppbar, npbar, +others?
 * */

double LowEnergySigma::BBbarTotal(int idA, int idB, double eCM) const {
  // @TODO: Reduce the total cross section if annihilation is impossible

  // Calculate effective energy, i.e. energy of protons with the same momenta
  double m0 = particleDataPtr->m0(2212);
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  double sBB = pow2(eCM);
  double sNN = 4 * pow2(m0) + (sBB - pow2(mA + mB)) 
                            * (sBB - pow2(mA - mB)) / sBB;
  double pLab = sqrt(sNN * (sNN - 4. * m0 * m0)) / (2. * m0);
  
  // Get parametrised cross section for ppbar
  double sigmaTotNN =
      (pLab < 0.3) ? 271.6 * exp(-1.1 * pLab * pLab)
    : (pLab < 5.)  ? 75.0 + 43.1 / pLab + 2.6 / pow2(pLab) - 3.9 * pLab
                   : HERAFit(38.4, 77.6, -0.64, 0.26, -1.2, pLab);

  // Scale by AQM factor
  // @TODO should this be elastic AQM factor?
  return sigmaTotNN * aqm(idA, idB) / aqmNN();
}

double LowEnergySigma::BBbarElastic(int idA, int idB, double eCM) const {
  // Calculate effective energy, i.e. energy of protons with the same momenta
  double m0 = particleDataPtr->m0(2212);
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  double sBB = pow2(eCM);
  double sNN = 4 * pow2(m0) + (sBB - pow2(mA + mB))
                            * (sBB - pow2(mA - mB)) / sBB;
  double pLab = sqrt(sNN * (sNN - 4. * m0 * m0)) / (2. * m0);

  // Get parametrised cross section for ppbar
  double sigmaElNN =
      (pLab < 0.3) ? 78.6
    : (pLab < 5.)  ? 31.6 + 18.3 / pLab - 1.1 / pow2(pLab) - 3.8 * pLab
                   : HERAFit(10.2, 52.7, -1.16, 0.125, -1.28, pLab);

  // Scale by AQM factor
  return sigmaElNN * aqm(idA, idB) / aqmNN();
}

double LowEnergySigma::BBbarNonDiff(int idA, int idB, double eCM) const {
  return BBbarTotal(idA, idB, eCM) - BBbarElastic(idA, idB, eCM) 
       - BBbarDiffractiveAX(idA, idB, eCM) - BBbarDiffractiveXB(idA, idB, eCM)
       - BBbarDiffractiveXX(idA, idB, eCM) - BBbarAnnihilation(idA, idB, eCM);
}

double LowEnergySigma::BBbarDiffractiveAX(int idA, int idB, double eCM) const {
  return BBDiffractiveAX(idA, -idB, eCM);
}

double LowEnergySigma::BBbarDiffractiveXB(int idA, int idB, double eCM) const {
  return BBDiffractiveXB(idA, -idB, eCM);
}

double LowEnergySigma::BBbarDiffractiveXX(int idA, int idB, double eCM) const {
  return BBDiffractiveXX(idA, -idB, eCM);
}

double LowEnergySigma::BBbarAnnihilation(int idA, int idB, double eCM) const {
  
  // Count number of each quark type
  vector<int> countA(3), countB(3);
  for (int quarksA = ( idA / 10) % 1000; quarksA > 0; quarksA /= 10)
    countA[quarksA % 10 - 1] += 1;
  for (int quarksB = (-idB / 10) % 1000; quarksB > 0; quarksB /= 10)
    countB[quarksB % 10 - 1] += 1;

  // Find the maximum number of simultaneous annihilations  
  int nMutual = 0;
  for (int i = 0; i < 3; ++i) 
    nMutual += min(countA[i], countB[i]);
  
  // Abort if no quarks can annihilate
  if (nMutual == 0)
    return 0.;

  // Calculate effective energy, i.e. energy of protons with the same momenta
  double m0 = particleDataPtr->m0(2212);
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  double sBB = pow2(eCM);
  double sNN = 4 * pow2(m0) + (sBB - pow2(mA + mB)) * (sBB - pow2(mA - mB)) / sBB;

  // Annihilation
  static constexpr double sigma0 = 120., A = 0.050, B = 0.6;
  double s0 = 4. * m0 * m0;
  double sigmaAnnNN = sigma0 * s0 / sNN
                    * ((A * A * s0) / (pow2(sNN - s0) + A * A * s0) + B);

  // Scale by AQM factor
  return sigmaAnnNN * aqm(idA, idB) / aqmNN();
}

//--------------------------------------------------------------------------

// Hadron-Meson section

/**TODO list for XM:
 * @TODO Diffraction using SaSDL
 * @TODO Consider parametrising ppiTotal, and define diff = total - res - elastic
 * @TODO Parametrise other cases explicitly, such as other pi+N and K+N
 * @TODO Get full elastic cross section (i.e. include resonant elastic)
 * @TODO Something with strangeness exchange?
 * @TODO UrQMD includes something about Danielewicz forward delay
 * @TODO Decide which resonances should be implemented
 * @TODO Check that all excited particles have the correct data, including id
 * @TODO Verify that mass-dependent widths are correct for all particles, esp. rho
 * @TODO Compare pi/K + p cross sections to PDG data
 * @TODO Compare a lot of cases to UrQMD
 **/

static Interpolator ppiStringData(1.9, 3.19642, {
    0., 0.597966, 1.6208, 2.64363, 3.66647, 4.6893, 5.71213, 6.67697,
    7.63029, 8.79865, 10.016, 11.3516, 12.4208, 13.3126, 13.984, 14.5885,
    15.0755, 15.4614, 15.7966, 16.0741, 16.3701, 16.6786, 16.9763,
    17.2395, 17.4756, 17.7065, 17.9797, 18.2733, 18.5514, 18.7887,
    18.9995, 19.1861, 19.3658, 19.4813, 19.5585, 19.6249, 19.6775,
    19.7497, 19.858, 20.0074, 20.1426, 20.2455, 20.3198, 20.3758,
    20.4241, 20.4542, 20.4502, 20.4389, 20.4559, 20.5093, 20.5774,
    20.6341
  });

static Interpolator ppiElData(1.975, 3.18545, 
  { 0., 1.75235, 2.18934, 2.44277, 2.63248, 2.78123, 2.92703, 3.10124,
    3.30385, 3.52057, 3.73107, 3.92351, 4.09365, 4.24155, 4.3692,
    4.47909, 4.57371, 4.65527, 4.72563, 4.78636, 4.83882, 4.9565,
    4.92306, 4.95654, 4.9852, 5.00957, 5.03016, 5.04739, 5.06161,
    5.07307, 5.08227, 5.08922, 5.09437, 5.09747, 5.09912, 5.09934,
    5.09824, 5.09596, 5.09174, 5.08824, 5.083, 5.07694, 5.07013, 5.06264,
    5.05453, 5.04584, 5.03664, 5.02696, 5.01684, 5.00633, 4.99546 }
);

double LowEnergySigma::XMTotal(int idX, int idM, double eCM) const {
  return XMNonDiffractive(idX, idM, eCM) + XMElastic(idX, idM, eCM) 
       + XMResonant(idX, idM, eCM);
}

double LowEnergySigma::XMNonDiffractive(int idX, int idM, double eCM) const {
  double sigmappi = ppiStringData(eCM);
  double aqmFactor = aqm(idX, idM) / aqm(2212, 211);
  return sigmappi * aqmFactor;
}

double LowEnergySigma::XMElastic(int idX, int idM, double eCM) const {
  if (particleDataPtr->isBaryon(idX)) 
    // Use parametrisation of ppi and scale by an aqm factor
    // @TBD: Should this scale by the elastic aqm or the total aqm?
    return ppiElData(eCM) * pow(aqm(idX, idM) / aqm(2212, 211), 2./3.);
  else
    // For meson-meson, return a fixed cross section
    return 5.;
}

double LowEnergySigma::XMResonant(int idX, int idM, double eCM) const {
  // For K_S and K_L, take average of K and Kbar
  if (idX == 310 || idX == 130)
    return 0.5 * (XMResonant(311, idM, eCM) + XMResonant(-311, idM, eCM));
  if (idM == 310 || idM == 130)
    return 0.5 * (XMResonant(idX, 311, eCM) + XMResonant(idX, -311, eCM));

  double sigmaRes = 0.;
  for (int idR : possibleResonances(idX, idM)) 
    sigmaRes += XMResonantPartial(idX, idM, idR, eCM);

  return sigmaRes;
}

int LowEnergySigma::getSignature(int baryonNumber, int charge, int nStrange) const {
  return 100 * baryonNumber
       +  10 * ((charge >= 0) ? charge : (10 + charge))
       +   1 * nStrange;
}

vector<int> LowEnergySigma::possibleResonances(int idX, int idM) const {
  ParticleDataEntry* entryA = particleDataPtr->findParticle(idX);
  ParticleDataEntry* entryB = particleDataPtr->findParticle(idM);
  int baryonNumber = entryA->isBaryon() + entryB->isBaryon();
  int charge = entryA->chargeType(idX) + entryB->chargeType(idM);
  int nStrange = entryA->nQuarksInCode(3) + entryB->nQuarksInCode(3);

  int signature = getSignature(baryonNumber, charge, nStrange);
  auto iter = signatureToParticles.find(signature);
  return (iter != signatureToParticles.end()) ? iter->second : vector<int>();
}

double LowEnergySigma::XMResonantPartial(int idX, int idM, int idR, 
    double eCM) const {

  double gammaR = hadronWidthsPtr->width(idR, eCM);
  double br = hadronWidthsPtr->branchingRatio(idR, idX, idM, eCM);
  
  if (gammaR == 0. || br == 0.)
    return 0.;

    // Find particle entries
  auto* entryR = particleDataPtr->findParticle(idR);
  auto* entryA = particleDataPtr->findParticle(idX);
  auto* entryB = particleDataPtr->findParticle(idM);

  // Calculate the resonance sigma
  double s = pow2(eCM), mA = entryA->m0(), mB = entryB->m0();
  double pCMS2 = 1 / (4 * s) * (s - pow2(mA + mB)) * (s - pow2(mA - mB));

  return GEVINVSQ2MB * M_PI / pCMS2
    * entryR->spinType() / (entryA->spinType() * entryB->spinType())
    * br * pow2(gammaR) / (pow2(entryR->m0() - eCM) + 0.25 * pow2(gammaR));
}

//==========================================================================

}