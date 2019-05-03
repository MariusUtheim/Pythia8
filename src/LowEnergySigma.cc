
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergySigma.h"

namespace Pythia8 {

// @TODO: Move this to a better place
static double clamp(double x, double min, double max) {
  return (x < min) ? min : (x > max) ? max : x;
}

//==========================================================================

// The LowEnergySigma class.

//--------------------------------------------------------------------------

// Returns int representing the overall process type:
//  0 - Collision not implemented (ids will not be ordered)
//  1 - BB
//  2 - BBbar
//  3 - XM
// 
// The canonical ordering of A and B satisfies two criteria:
//   1) |A| >= |B|, and 2) A > 0.
// 
// Implications:
//  - If both A and B are negative, they are replaced by their antiparticles
//  - In BBbar, A is the particle and B is the antiparticle
//  - In XM, B is always a meson and X is a baryon (not antibaryon) or meson 
//
int LowEnergySigma::canonicalForm(int& idA, int& idB) const {
  if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB))
    return 0;

  // Ensure |A| >= |B|
  if (abs(idA) < abs(idB))
    swap(idA, idB);

  // Ensure A > 0
  if (idA < 0) { 
    idA = -idA;
    if (particleDataPtr->hasAnti(idB))
      idB = -idB; 
  }

  // Get id of overall collision type
  if (particleDataPtr->isMeson(idB))
    return 3; // XM
  else if (idB < 0)
    return 2; // BBbar
  else
    return 1; // BB
}

//--------------------------------------------------------------------------

double LowEnergySigma::sigmaTotal(int idA, int idB, double eCM) const {
  if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB))
    return 0.; // @TODO Error message in this case?

  // If energy is less than the hadron masses, return 0.
  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return 0.;

  switch (canonicalForm(idA, idB)) {
    case 0: return 0.; // @TODO: Probably give an error message?
    case 1: return sigmaTotalBB(idA, idB, eCM);
    case 2: return sigmaTotalBBbar(idA, idB, eCM);
    case 3: return sigmaTotalXM(idA, idB, eCM);
    default: return 0.; // @TODO: This should never occur
  }
}

//--------------------------------------------------------------------------

map<int, double> LowEnergySigma::sigmaPartial(int idA, int idB, double eCM) const {
  if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB))
    return map<int, double>(); // @TODO Error message in this case?

  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return map<int, double>();

  switch (canonicalForm(idA, idB)) {
    case 1: {
      double tot = sigmaTotalBB(idA, idB, eCM);
      double el  = sigmaElasticBB(idA, idB, eCM);
      return map<int, double>{
        { 0, tot },
        { 1, tot - el },
        { 2, el }
      };
    }

    case 2:
      return sigmaPartialBBbar(idA, idB, eCM);
    
    case 3: {
      map<int, double> result;
      double res = 0.;
      for (auto idR : lowEnergyResPtr->getPossibleResonances(idA, idB)) {
        double partial = lowEnergyResPtr->getPartialResonanceSigma(idA, idB, idR, eCM);
        if (partial > 0) {
          res += partial;
          result.emplace(idR, partial);
        }
      }
      double inel = sigmaInelXM(idA, idB, eCM);
      double elastic = sigmaElasticXM(idA, idB, eCM);
      result.emplace(0, res + inel + elastic);
      result.emplace(1, inel);
      result.emplace(2, elastic);
      result.emplace(7, res);
      return result;
    }

    default: // @TODO: This should not be possible to reach
      return map<int, double>(); 
  }
}

//--------------------------------------------------------------------------

double LowEnergySigma::sigmaPartial(int idA, int idB, double eCM, int proc) const {
  // Note: A shorthand way of calculating this would be to get the map of
  // partial cross sections, then return the value if the specified process is
  // contained in the map. This risks being slow however, e.g. when calculating
  // the cross section for a particular resonance, we don't need to calculate
  // the cross sections for all other resonances.

  if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB))
    return 0.; // @TODO Error message in this case?

  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return 0.;

  switch (canonicalForm(idA, idB)) {
    case 0: return 0.; // @TODO: Probably give an error message?

    case 1: // BB
      sigmaSaSDL.calcDiff(idA, idB, eCM * eCM,
                          particleDataPtr->m0(idA), particleDataPtr->m0(idB));
      switch (proc) {

        case 0: // Inelastic
          return sigmaTotalBB(idA, idB, eCM) - sigmaElasticBB(idA, idB, eCM);

        case 2: 
          //sigmaSaSDL.calcTotEl(idA, idB, eCM * eCM,
          //                particleDataPtr->m0(idA), particleDataPtr->m0(idB));
          //return sigmaSaSDL.sigEl;
          return sigmaElasticBB(idA, idB, eCM);

        case 1: case 3: case 4: case 5: {
          sigmaSaSDL.calcTotEl(idA, idB, eCM * eCM,
                         particleDataPtr->m0(idA), particleDataPtr->m0(idB));
          sigmaSaSDL.calcDiff(idA, idB, eCM * eCM,
                         particleDataPtr->m0(idA), particleDataPtr->m0(idB));
          double sigND = sigmaSaSDL.sigTot - sigmaSaSDL.sigEl
            - sigmaSaSDL.sigAX - sigmaSaSDL.sigXB - sigmaSaSDL.sigAXB - sigmaSaSDL.sigXX;
          double factor = proc == 1 ? sigND
                        : proc == 3 ? sigmaSaSDL.sigAX 
                        : proc == 4 ? sigmaSaSDL.sigXB
                                    : sigmaSaSDL.sigXX;
          factor /= (sigmaSaSDL.sigTot - sigmaSaSDL.sigEl);
        
          return factor * (sigmaTotalBB(idA, idB, eCM) - sigmaElasticBB(idA, idB, eCM));
        }

        default: return 0;
      }

    case 2: { // BBbar
      auto sigmas = sigmaPartialBBbar(idA, idB, eCM);
      auto iter = sigmas.find(proc);
      return (iter == sigmas.end()) ? 0. : iter->second;
    }

    case 3: 
      switch (proc) {
        case 1: return sigmaInelXM(idA, idB, eCM);
        case 2: return sigmaElasticXM(idA, idB, eCM);
        case 7: return lowEnergyResPtr->getResonanceSigma(idA, idB, eCM);
        default:
          return (abs(proc) <= 100) ? 0. 
              : lowEnergyResPtr->getPartialResonanceSigma(idA, idB, proc, eCM);
      }

    default: // @TODO: This should not be possible to reach
      return 0; 
  }
}

//--------------------------------------------------------------------------

double LowEnergySigma::aqm(int idA, int idB) const {
  double mesA = particleDataPtr->isMeson(idA);
  double mesB = particleDataPtr->isMeson(idB);
  return 40 * pow(2./3., mesA + mesB)
       * (1 - 0.4 * abs(particleDataPtr->nStrangeQuarks(idA)) / (mesA ? 2 : 3))
       * (1 - 0.4 * abs(particleDataPtr->nStrangeQuarks(idB)) / (mesB ? 2 : 3));
}

double LowEnergySigma::aqmNN() const {
  return 40;
}

static double ReggeFit(double z, double y1, double y2, double s) {
  return z + 0.308 * pow2(log(s / 28.998)) + y1 * pow(s, -0.458) - y2 * pow(s, -0.545);
}

static double HERAFit(double a, double b, double n, double c, double d, double p) {
  return a + b * pow(p, n) + c * pow2(log(p)) + d * log(p);
}

//--------------------------------------------------------------------------

// Baryon-Baryon section

/**@TODO list for BB:
 *  Check that the tables are correct and sufficiently smooth, compare w/ data
 *  Implement strangeness exchange
 *  Implement parametrisation for Lambda+p and Sigma+p special cases
 *  Implement D+N and D+D collisions
 *  Decide when diffraction occurs and when strings are formed
 *  Do something abour charm and bottom? 
 **/

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

// Strangeness exchange
double LowEnergySigma::sigmaStrEx(int, int, double) const {
  // @TODO Implement this
  return 0.;
}

double LowEnergySigma::sigmaTotalBB(int idA, int idB, double eCM) const {
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
    return aqm(idA, idB) + sigmaStrEx(idA, idB, eCM);
  }
}

double LowEnergySigma::sigmaElasticBB(int idA, int idB, double eCM) const {
  // Fit pp/nn/pn
  if ((idA == 2112 || idA == 2212) && (idB == 2112 || idB == 2212)) {
    if (eCM < 2 * particleDataPtr->m0(2212) + particleDataPtr->m0(211))
      return sigmaTotalBB(idA, idB, eCM);

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
    // @TODO: UrQMD puts a threshold on this value, returning 0 if < 1 MeV
    return 0.039 * pow(aqm(idA, idB), 2./3.);
  }
}

//--------------------------------------------------------------------------

// Baryon-Antibaryon section

/**@TODO list for BBbar:
 *  Check that sNN is correct
 *  UrQMD actually uses Regge fit instead of HERA fit for sigmaTotNN
 *  sigmaTotNN and sigmaElNN do not match data well for pLab < 0.3
 *  Should there be a different parametrisation for npbar?
 *  Check that the aqmFactor is correct
 *  Figure out if the annihilation cross section parametrisation is up to date
 *  Decide how to split the diffractive cross sections among specific processes
 *  Verify that all the process numbers are actually correct
 * Compare with data cases: ppbar, npbar, +others?
 * */

double LowEnergySigma::sigmaTotalBBbar(int idA, int idB, double eCM) const {
  
  // Calculate effective energy, i.e. energy of protons with the same momenta
  double m0 = particleDataPtr->m0(2212);
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  double sBB = pow2(eCM);
  double sNN = 4 * pow2(m0) + (sBB - pow2(mA + mB)) * (sBB - pow2(mA - mB)) / sBB;
  double pLab = sqrt(sNN * (sNN - 4. * m0 * m0)) / (2. * m0);

  // Get parametrised cross section for ppbar
  double sigmaTotNN =
      (pLab < 0.3) ? 271.6 * exp(-1.1 * pLab * pLab)
    : (pLab < 5.)  ? 75.0 + 43.1 / pLab + 2.6 / pow2(pLab) - 3.9 * pLab
                   : HERAFit(38.4, 77.6, -0.64, 0.26, -1.2, pLab);

  // Get scale factor (from AQM)
  double aqmFactor = aqm(idA, idB) / aqmNN();
  
  return sigmaTotNN * aqmFactor;
}

// @TODO: Consider making separate functions. For now I chose to put them all
// in one function, since many partial cross sections depend on the same initial
// computations, such as sNN and sigmaTotal
map<int, double> LowEnergySigma::sigmaPartialBBbar(int idA, int idB, double eCM) const {
  // Calculate effective energy, i.e. energy of protons with the same momenta
  double m0 = particleDataPtr->m0(2212);
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);
  double sBB = pow2(eCM);
  double sNN = 4 * pow2(m0) + (sBB - pow2(mA + mB)) * (sBB - pow2(mA - mB)) / sBB;
  double pLab = sqrt(sNN * (sNN - 4. * m0 * m0)) / (2. * m0);

  // Get parametrised cross section for ppbar
  double sigmaTotNN =
      (pLab < 0.3) ? 271.6 * exp(-1.1 * pLab * pLab)
    : (pLab < 5.)  ? 75.0 + 43.1 / pLab + 2.6 / pow2(pLab) - 3.9 * pLab
                   : HERAFit(38.4, 77.6, -0.64, 0.26, -1.2, pLab);

  // Elastic
  double sigmaElNN =
      (pLab < 0.3) ? 78.6
    : (pLab < 5.)  ? 31.6 + 18.3 / pLab - 1.1 / pow2(pLab) - 3.8 * pLab
                   : HERAFit(10.2, 52.7, -1.16, 0.125, -1.28, pLab);
  
  // Annihilation
  double sigma0 = 120., A = 0.050, B = 0.6;
  double s0 = 4. * m0 * m0;
  double sigmaAnnNN = sigma0 * s0 / sNN
                    * ((A * A * s0) / (pow2(sNN - s0) + A * A * s0) + B);

  // Diffractive (string + inelastic)
  double sigmaInelasticNN = sigmaTotNN - sigmaElNN - sigmaAnnNN;
  double t = clamp((eCM - 3.) / (5. - 3.), 0., 1.);
  double sigmaStringNN = t * sigmaInelasticNN;
  double sigmaDiffNN = (1 - t) * sigmaInelasticNN;

  // Get scale factor (from AQM)
  double aqmFactor = aqm(idA, idB) / aqmNN();

  return map<int, double>{
    { 0, sigmaTotNN    * aqmFactor },
    { 1, sigmaStringNN * aqmFactor },
    { 2, sigmaElNN     * aqmFactor },
    { 3, sigmaDiffNN   * aqmFactor },
    // @TODO: Also other diffractive cases
    { 6, sigmaAnnNN    * aqmFactor }
  };
}

//--------------------------------------------------------------------------

// Hadron-Meson section

/**@TODO list for XM:
 *  Consider parametrising ppiTotal, and define diff = total - res - elastic
 *  Parametrise other cases explicitly, such as other pi+N and K+N
 *  Get full elastic cross section (i.e. include resonant elastic)
 *  Something with strangeness exchange?
 *  UrQMD includes something about Danielewicz forward delay
 *  Decide which resonances should be implemented
 *  Check that all excited particles have the correct data, including id
 *  Verify that mass-dependent widths are correct for all particles
 *  Compare pi/K + p cross sections to PDG data
 *  Compare a lot of cases to UrQMD
 **/

static Interpolator ppiDiffData(1.9, 3.19642, {
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
    5.05453, 5.04584, 5.03664, 5.02696, 5.01684, 5.00633, 4.99546 });


// Total = resonant + elastic + diffractive(including strings)
double LowEnergySigma::sigmaTotalXM(int idX, int idM, double eCM) const {
  return lowEnergyResPtr->getResonanceSigma(idX, idM, eCM)
       + sigmaElasticXM(idX, idM, eCM) + sigmaInelXM(idX, idM, eCM);
}

double LowEnergySigma::sigmaInelXM(int idX, int idM, double eCM) const {
  double sigmaDiffppi = ppiDiffData(eCM);
  double aqmFactor = aqm(idX, idM) / aqm(2212, 211);
  return sigmaDiffppi * aqmFactor;
}

double LowEnergySigma::sigmaElasticXM(int idX, int idM, double eCM) const {
  if (particleDataPtr->isBaryon(idX))
    // Return parametrisation of ppi, scaled by an aqm factor
    // @TODO: Should this scale by the elastic aqm or the total aqm?
    return ppiElData(eCM) * pow(aqm(idX, idM) / aqm(2212, 211), 2./3.);
  else
    return 5.;
}

//==========================================================================

}