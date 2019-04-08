
#include "Pythia8/Interpolator.h"
#include "Pythia8/LowEnergySigma.h"

namespace Pythia8 {

double LowEnergySigma::sigmaTotal(int idA, int idB, double eCM) const {
  if (!particleDataPtr->isHadron(idA) || !particleDataPtr->isHadron(idB))
    return 0.; // @TODO Error message in this case?

  if (particleDataPtr->m0(idA) + particleDataPtr->m0(idB) > eCM)
    return 0.; 

  int bA = particleDataPtr->baryonNumberType(idA), 
      bB = particleDataPtr->baryonNumberType(idB);

  // If two baryons
  if (bA && bB) {
    // If two antibaryons, negate both
    if (bA < 0 && bB < 0)
    { idA = -idA; idB = -idB; }

    // If one antibaryon
    if (bA * bB < 0) {
      // Ensure idA is the baryon and idB is antibaryon
      if (bA < 0)
        swap(idA, idB);

      // Get BBbar total cross section
      return sigmaTotalBBbar(idA, idB, eCM);
    }
    // If no antibaryons (or both are antibaryons and ids have been negated)
    else {
      // Get BB total cross section
      return sigmaTotalBB(idA,  idB, eCM);
    }
  }
  else {
    // If there is an antibaryon, negate both
    // @TODO Check if both particles have an antiparticle?
    if (bA < 0 || bB < 0)
    { idA = -idA; idB = -idB; }

    // If there is a baryon, ensure it is A
    if (bB != 0)
      swap(idA, idB);

    // Get XM total cross section
    return sigmaTotalXM(idA, idB, eCM);
  }
}

double LowEnergySigma::sigmaPartial(int idA, int idB, double eCM, int proc) const {
  switch (proc) {
    case 7:
      return lowEnergyResPtr->getResonanceSigma(idA, idB, eCM);

    default:
      if (particleDataPtr->isHadron(proc))
        return lowEnergyResPtr->getPartialResonanceSigma(idA, idB, proc, eCM);
  }

  return 0.;
}

static Interpolator ppTotalData(1.88, 5.0, {
  335.561, 99.5353, 32.9358, 27.753, 24.4147, 23.7205, 23.8078,
  24.1367, 24.0177, 25.2023, 26.1183, 29.1103, 31.7825, 34.3526, 
  38.1058, 41.6417, 43.6273, 45.439, 46.5519, 47.5467, 47.5348, 
  47.5201, 47.4971, 47.4584, 47.3441, 47.2299, 47.2413, 47.2926, 
  47.2723, 47.1574, 47.0425, 46.9276, 46.7771, 46.5213, 46.2655, 
  46.0097, 45.7539, 45.4981, 45.2562, 45.0241, 44.792, 44.5599, 
  44.3278, 44.3162, 44.3281, 44.34, 44.2698, 44.1811, 44.0924, 43.9772,
  43.7315, 43.4857, 43.24, 42.9942, 42.7736, 42.6152, 42.4568, 42.2983, 
  42.1399, 41.9815, 41.8231, 41.7054, 41.6882, 41.6711, 41.6539, 
  41.6367, 41.6195, 41.6023, 41.5851, 41.5679, 41.5507, 41.5335, 
  41.5164, 41.4992, 41.482, 41.4648, 41.4476, 41.4304, 41.4132, 41.396,
  41.3788, 41.3617, 41.3445, 41.3273, 41.3101, 41.2929, 41.2757,
  41.2585, 41.2413, 41.2241, 41.207, 41.1898, 41.1726, 41.1554, 
  41.1382, 41.121, 41.1038, 41.0866, 41.0694, 41.0522, 41.0351, 
  41.0179, 41.0007, 40.9835, 40.9663, 40.9491, 40.9319, 40.9147,
  40.8975, 40.8804, 40.8632, 40.846, 40.8288, 40.8116, 40.7944, 
  40.7772, 40.76, 40.7428, 40.7257, 40.7085, 40.6913, 40.6741, 40.6569,
  40.6397, 40.6225, 40.6053, 40.5881, 40.571, 40.5538, 40.5366, 
  40.5194, 40.5022, 40.485, 40.4678, 40.4506, 40.4334, 40.4162, 
  40.3991, 40.3819, 40.3647, 40.3475, 40.3303, 40.3131, 40.2959,
  40.2787, 40.2615, 40.2444, 40.2272, 40.21, 40.1928, 40.1756, 40.1584,
  40.1412, 40.124, 40.1068, 40.0897, 40.0725
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

static double ReggeFit(double z, double y1, double y2, double s) {
  return z + 0.308 * pow2(log(s / 28.998)) + y1 * pow(s, -0.458) - y2 * pow(s, -0.545);
}

double LowEnergySigma::sigmaTotalBB(int idA, int idB, double eCM) const {
// pp tot: z = 35.45, y1 = 42.53, y2 = 33.34 
// pn tot: z = 35.80, y1 = 40.15, y2 = 30.00
// pp/pn el: (use HERA)

  // Look for parametrisation
  if ((idA == 2212 && idB == 2212)
   || (idA == 2112 && idB == 2112)) {
    double t = (eCM - 3.) / (5. - 3.);
    double parametrised = ReggeFit(35.45, 42.53, 33.34, eCM * eCM);
    return (1 - t) * ppTotalData(eCM) + t * parametrised;
  }
  if (idA == 2212 && idB == 2112)
  {
    double t = (eCM - 3.) / (5. - 3.);
    double parametrised = ReggeFit(35.80, 40.15, 30.00, eCM * eCM);
    return (1 - t) * pnTotalData(eCM) + t * parametrised;

  }
  // @TODO: Something special for Delta1232+N or Delta1232+Delta1232
  else {
    double sigmaAQM = 40. * (1 - 0.4 * strangenessFactor(idA)) 
                          * (1 - 0.4 * strangenessFactor(idB));

    // @TODO Add strangeness exchange for Lambda+Sigma or Xi+N
    {
      double sigmaSTREX = 0.; // @TODO get STREX
      sigmaAQM += sigmaSTREX; 
    }

    return sigmaAQM;
  }
  

  return 0.;
}

double LowEnergySigma::sigmaTotalBBbar(int idA, int idB, double eCM) const {
  // Assumption: idA > 0, idB < 0

  double m0 = particleDataPtr->m0(2212);
  double mA = particleDataPtr->m0(idA), mB = particleDataPtr->m0(idB);

  // Calculate effective energy, i.e. energy of protons with the same momenta
  double sBB = pow2(eCM);
  double sNN = pow2(m0) + (sBB - pow2(mA + mB)) * (sBB - pow2(mA - mB)) / sBB;
  double pLab = sqrt(sNN * (sNN - 4. * m0 * m0)) / (2. * m0);

  // Get parametrised cross section for ppbar
  double sigmaTotNN =
      (pLab < 0.3) ? 271.6 * exp(-1.1 * pLab * pLab)
    : (pLab < 5.)  ? 75.0 + 43.1 / pLab + 2.6 / pow2(pLab) - 3.9 * pLab
                   : 38.4 + 77.6 * pow(pLab, -0.64) 
                          + 0.26 * pow2(log(pLab)) - 1.2 * log(pLab);

  // Get scale factor (from AQM)
  double aqmFactor;
  if ((idA == 2212 && idB == -2212) || (idA == 2112 && idB == -2112)
   || (idA == 2212 && idB == -2112) || (idA == 2112 && idB == -2212))
    aqmFactor = 1.;
  else
    aqmFactor = (1 - 0.4 / 3 * particleDataPtr->strangeness(idA))
              * (1 - 0.4 / 3 * particleDataPtr->strangeness(idB));
  
  return sigmaTotNN * aqmFactor;

  //// PARTIAL CROSS SECTIONS

//  // Elastic
//  double sigmaElNN =
//      (pLab < 0.3) ? 78.6
//    : (pLab < 5.)  ? 31.6 + 18.3 / pLab - 1.1 / pow2(pLab) - 3.8 * pLab
//                   : 10.2 + 52.7 * pow(pLab, -1.16) 
//                          + 0.125 * pow2(log(pLab)) - 1.28 * log(pLab);
//  
//  // Annihilation
//  double sigma0 = 120., Asq = 7.7404804e-5, B = 0.6;
//  double s0 = 4. * m0 * m0;
//  double sigmaAnnNN = sigma0 * s0 / sNN
//                    * (Aqs / (pow2(sNN - s0) + Asq) + B);
//
//  // Diffractive (string + inelastic)
//  double sigmaDiffNN = sigmaTotNN - sigmaElNN - sigmaAnnNN;
//  double t = median(0, 1, (eCM - 3.) / (5. - 3.));
//  double sigmaStringNN = t * sigmaDiffNN;
//  double sigmaInelasticNN = (1 - t) * sigmaDiffNN;
//
}

//vector<pair<double, double>> ppiDiffData = 
//  {{1.8, 0.}, {1.9, 0.}, {2.01298, 4.08957}, {2.05891, 5.93332},
//  {2.10384, 7.60645}, {2.14785, 9.62448}, {2.19099, 11.8857}, {2.2333, 
//  13.3665}, {2.27484, 14.4379}, {2.31563, 15.2175}, {2.35573, 15.7593}, 
//  {2.39516, 16.1886}, {2.43395, 16.6593}, {2.47214, 17.1056}, {2.50975, 
//  17.4548}, {2.54681, 17.7907}, {2.58334, 18.2172}, {2.61936, 18.6103}, 
//  {2.6549, 18.9261}, {2.68997, 19.1831}, {2.72459, 19.4272}, {2.75877, 
//  19.5375}, {2.79254, 19.6255}, {2.82591, 19.6942}, {2.85889, 19.7995}, 
//  {2.89149, 19.9928}, {2.92374, 20.1639}, {2.95563, 20.2851}, {2.98718, 
//  20.3579}, {3.0184, 20.4203}, {3.0493, 20.4568}, {3.07989, 20.4481}, 
//  {3.11019, 20.4328}, {3.14019, 20.4847}, {3.16991, 20.5672}, {3.19642, 
//  20.6261}};

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

double LowEnergySigma::sigmaTotalXM(int idX, int idM, double eCM) const {
  double sigmaRes = lowEnergyResPtr->getResonanceSigma(idX, idM, eCM);

  // @TODO Something special for K_S and K_L

  // Get parametrised Npi cross section
  double sigmaElppi = ppiElData(eCM), sigmaDiffppi = ppiDiffData(eCM);

  if ((idX == 2212 || idX == 2112)
      && (idM == 211 || idM == 111 || idM == -211))
    return sigmaRes + sigmaElppi + sigmaDiffppi;
  
  // Total non-resonance cross section
  double sigmaNonresppi = sigmaElppi + sigmaDiffppi;
  
  // Correction factor from AQM
  // @TODO Check that this is correct
  double AQMfactor = 
      (particleDataPtr->isMeson(idX) ? 2./3. : 1.)
    * strangenessFactor(idX) * strangenessFactor(idM);
  
  return sigmaRes + AQMfactor * sigmaNonresppi;
}

}