// New parametrizations for sigma_diffractive
// by Appleby, Barlow, Molson, Serluca and Toader (ABMST).

//==========================================================================

#include "Pythia8/Basics.h"
using namespace Pythia8;

// Common values.
const double SPROTON    = 0.8803544;
const double MPROTON    = 0.9382720;
const double MPINEU     = 0.1349766;

//==========================================================================

// Single diffractive cross section according to
// Appleby, Barlow, Molson, Serluca and Toader (ABMST).

const double SDAI[] = { 0.624529, 3.09088, 4.0, 177.217};
const double SDBI[] = { 2.5835, 4.51487, 3.03392, 5.86474};
const double SDCI[] = { 0.0, 0.186211, 10.0, 21.0029};
const double EPSI[] = { 0.08, -0.4525};
const double ALPP[] = { 0.25, 0.93};
const double POMS[] = { -0.25, -1.15, -0.05, 0.4, 0.5, 0.4597, 5.7575};
const double PIPS[] = { 13.63, 0.0808, 31.79, -0.4525, 14.4 };
const double MRES[] = { 1.44, 1.52, 1.68, 2.19};
const double WRES[] = { 0.325, 0.130, 0.140, 0.450};
const double CRES[] = { 3.07, 0.4149, 1.108, 0.9515};
const double FORM[] = { 2.79, 1., 0.71};

//--------------------------------------------------------------------------

// dsigma_sd/(dt dxi)(s, t, xi).

double dsigSDdtdxi( double s, double t, double xi, int iPart = 0) {

  // Vanishing cross section below the p + pi threshold.
  double spipSum = pow2( MPROTON + MPINEU);
  double spipDif = pow2( MPROTON - MPINEU);
  double m2X     = xi * s;
  if (m2X < spipSum) return 0.;

  // Calculate t range. Vanishing cross section if outside that range.
  double mu1     = SPROTON / s;
  double mu3     = xi;
  double rootv   = (1. - 4. * mu1) * (pow2(1. - mu1 - mu3) - 4. * mu1 * mu3);
  if (rootv <= 0.) return 0.;
  double tMin    = -0.5 * s * ( 1. - 3. * mu1 - mu3 + sqrt(rootv) );
  double tMax    = s * s * mu1 * pow2(mu3 - mu1) / tMin;
  if (t < tMin || t > tMax) return 0.;

  // Separation between low- and high-mass diffraction.
  // For low-mass evaluate P+R terms at the cut scale to allow matching.
  double m2Cut   = pow2( (s < 4000.) ? 3. : 3. + 0.6 * log( s / 4000.) );
  bool isHighM   = (m2X > m2Cut);
  double xiHigh  = (isHighM) ? xi : m2Cut / s;

  // Value of Pomeron/Reggeon trajectory.
  double alp0Pom = 1. + EPSI[0];
  double alp0Reg = 1. + EPSI[1];
  double alptPom = alp0Pom + ALPP[0] * t;
  double alptReg = alp0Reg + ALPP[1] * t;

  // PPP term, split by t range.
  /*
  double dsigPPP = t / (t + POMS[2])
    * pow( xiHigh, alp0Pom - 2. * alptPom) * pow( s, EPSI[0]);
  if (t > POMS[0]) dsigPPP *= POMS[3] + POMS[4] * t;
  else             dsigPPP *= (SDAI[0] * exp(SDBI[0] * t) + SDCI[0]);
  if (t < POMS[1]) dsigPPP *=  1. + POMS[5] * (POMS[1] - t)
    + POMS[6] * pow2(POMS[1] - t);
  */
  double dsigPPP = pow( xiHigh, alp0Pom - 2. * alptPom) * pow( s, EPSI[0]);
  if (t > POMS[0]) dsigPPP *= POMS[3] + POMS[4] * t;
  else dsigPPP *= ((SDAI[0] * exp(SDBI[0] * t) + SDCI[0])) * t / (t + POMS[2]);
  if (t < POMS[1]) dsigPPP *=  1. + POMS[5] * (POMS[1] - t)
    + POMS[6] * pow2(POMS[1] - t);

  // PPR, RRP and RRR terms.
  double dsigPPR = (SDAI[1] * exp(SDBI[1] * t) + SDCI[1])
    * pow( xiHigh, alp0Reg - 2. * alptPom) * pow( s, EPSI[1]);
  double dsigRRP = (SDAI[2] * exp(SDBI[2] * t) + SDCI[2])
    * pow( xiHigh, alp0Pom - 2. * alptReg) * pow( s, EPSI[0]);
  double dsigRRR = (SDAI[3] * exp(SDBI[3] * t) + SDCI[3])
    * pow( xiHigh, alp0Reg - 2. * alptReg) * pow( s, EPSI[1]);

  // Pion exchange term.
  double fForm   = (4. * SPROTON - FORM[0] * t)
                 / ( (4. * SPROTON - FORM[1] * t) * pow2(1. - t / FORM[2]) );
  double alptPi  = ALPP[1] * (t - pow2(MPINEU));
  double sigpip  = PIPS[0] * pow( xiHigh * s, PIPS[1])
                 + PIPS[2] * pow( xiHigh * s, PIPS[3]);
  double dsigPi  = PIPS[4] / (4. * M_PI) * (-t) / pow2(t - pow2(MPINEU))
                 * pow2(fForm) * pow( xiHigh, 1. - 2. * alptPi) * sigpip;

  // Done for high-mass diffraction.
  double dsigBkg = 0.;
  double dsigRes = 0.;
  double dsigTot = 0.;
  if (isHighM) dsigTot = dsigPPP + dsigPPR + dsigRRP + dsigRRR + dsigPi;

  // Low-mass diffraction: smoothly dampen high-mass contribution.
  else {
    double dsigCut = dsigPPP + dsigPPR + dsigRRP + dsigRRR + dsigPi;
    double derdCut = ( (alp0Pom - 2. * alptPom) * dsigPPP
      + (alp0Reg - 2. * alptPom) * dsigPPR
      + (alp0Pom - 2. * alptReg) * dsigRRP
      + (alp0Reg - 2. * alptReg) * dsigRRR
      + (1. - 2. * alptPi) * dsigPi ) / xiHigh;
    double xiThr   = spipSum / s;
    double derdDif = derdCut * (xiHigh - xiThr) - dsigCut;
    double coef2   = derdDif / pow2(xiHigh - xiThr);
    double coef1   = derdCut - 2. * derdDif / (xiHigh - xiThr);
    dsigBkg        = coef2 * pow2(xi - xiThr) + coef1 * (xi - xiThr);

    // Add low-mass resonances.
    double dsigSub = 0.;
    double qRef    = 0.5 * sqrt( (m2X - spipSum) * (m2X - spipDif) / m2X);
    for (int i = 0; i < 4; ++i) {
      double m2Now = pow2( MRES[i]);
      double qNow  = 0.5 * sqrt( (m2Now - spipSum) * (m2Now - spipDif) / m2Now);
      double mwNow = MRES[i] * WRES[i] * pow( qRef / qNow, 2 * i + 3)
                   * pow( (1. + 5. * qNow) / (1. + 5. * qRef), i + 1);
      dsigRes     += CRES[i] * mwNow / (pow2(m2X - m2Now) + pow2(mwNow));
      dsigSub     += CRES[i] * mwNow / (pow2(m2Cut - m2Now) + pow2(mwNow));
   }
    // Include Jacobian and t dependence, and subtract to vanish at cut.
    dsigSub       *= (xi - xiThr) / (xiHigh - xiThr);
    dsigRes        = (dsigRes / xi - dsigSub / xiHigh) * exp(13.5 * (t + 0.05));

    // Total low-mass contribution.
    dsigTot        = dsigBkg + dsigRes;
  }

  // Return sum, or for debug individual terms.
  if      (iPart == 0) return dsigTot;
  else if (iPart == 1 && isHighM) return dsigPPP;
  else if (iPart == 2 && isHighM) return dsigPPR;
  else if (iPart == 3 && isHighM) return dsigRRP;
  else if (iPart == 4 && isHighM) return dsigRRR;
  else if (iPart == 5 && isHighM) return dsigPi;
  else if (iPart == 6) return dsigBkg;
  else if (iPart == 7) return dsigRes;
  else                 return 0.;
}

//--------------------------------------------------------------------------

// Integrate over dt to obtain dsigma_sd/dxi(s, xi) in range tMin < t < tMax.

double dsigSDdxi( double s, double xi, double tMinIn = -1e10,
  double tMaxIn = 0., int iPart = 0) {

  // Number of integration points.
  int nPoints  = 200;

  // Calculate t range.
  double mu1   = SPROTON / s;
  double mu3   = xi;
  double rootv = (1. - 4. * mu1) * (pow2(1. - mu1 - mu3) - 4. * mu1 * mu3);
  if (rootv <= 0.) return 0.;
  double tMin  = -0.5 * s * ( 1. - 3. * mu1 - mu3 + sqrt(rootv) );
  double tMax  = s * s * mu1 * pow2(mu3 - mu1) / tMin;

  // Impose further t constraints from input. Check if range closed.
  tMin = max( tMin, tMinIn);
  tMax = min( tMax, tMaxIn);
  if (tMin > tMax) return 0.;

  // Prepare integration.
  double slope = -0.5 * log(xi);
  double etMin = exp(slope * tMin);
  double etMax = exp(slope * tMax);

  // Do integration by uniform steps in exp(slope * t).
  double dsigdxi= 0.;
  double y, yt, t;
  for (int i = 0; i < nPoints; ++i) {
    y          = (i + 0.5) / nPoints;
    yt         = etMin + y * (etMax - etMin);
    t          = log(yt) / slope;
    dsigdxi   += dsigSDdtdxi(s, t, xi, iPart) / yt;
  }

  // Normalize and done.
  dsigdxi *= (etMax - etMin) / (nPoints * slope);
  return dsigdxi;

}

//--------------------------------------------------------------------------

// Integrate over xi to obtain dsigma_sd/dt(s) in range xi_min < xi < xi_max.

double dsigSDdt( double s, double t, double xiMinIn = 0., double xiMaxIn = 1.,
  int iPart = 0) {

  // Density of integration points in ln(xi) for xi < xiDiv or else in xi.
  double xiDiv    = 0.1;
  double dlnxiRaw = 0.1;
  double dxiRaw   = 0.01;
  double mXmin    = MPROTON + MPINEU;

  // Approximate restrictions on xi range. Check it is not closed.
  double dsigdt   = 0.;
  double xiMin    = max( xiMinIn, mXmin * mXmin / s);
  double xiMax    = min( xiMaxIn, sqrt( -t / SPROTON) );
  if (xiMin >= xiMax) return 0.;

  // Integration in xi: check size of affected range and adjust dxi.
  if (xiMax > xiDiv) {
    double xiMinP = max( xiDiv, xiMin);
    int    nxi    = 2 + (xiMax - xiMinP) / dxiRaw;
    double dxi    = (xiMax - xiMinP) / nxi;
    for (int ixi = 0; ixi < nxi; ++ixi)
      dsigdt += dxi * dsigSDdtdxi( s, t, xiMinP + dxi * (ixi + 0.5), iPart);
  }

  // Integration in ln(xi): check size of affected range and adjust dlnxi.
  if (xiMin < xiDiv) {
    double xiMaxP = min( xiDiv, xiMax);
    int    nlnxi  = 2 + log( xiMaxP / xiMin) / dlnxiRaw;
    double dlnxi  = log( xiMaxP / xiMin) / nlnxi;
    for (int ilnxi = 0; ilnxi < nlnxi; ++ilnxi) {
      double xi   = xiMin * exp( dlnxi * (ilnxi + 0.5));
      dsigdt     += dlnxi * xi * dsigSDdtdxi( s, t, xi, iPart);
    }
  }

  // Done.
  return dsigdt;
}

//--------------------------------------------------------------------------

// Integrate over t and xi to obtain sigma_sd(s) in range  tMin < t < tMax
// and xi_min < xi < xi_max.

double sigSD( double s, double tMinIn = -1e10, double tMaxIn = 0.,
  double xiMinIn = 0., double xiMaxIn = 1., int iPart = 0) {

  // Density of integration points in ln(xi) for xi < xiDiv or else in xi.
  double xiDiv    = 0.1;
  double dlnxiRaw = 0.1;
  double dxiRaw   = 0.01;
  double mXmin    = 1.2;

  // Restrictions on xi range. Check it is not closed.
  double sig      = 0.;
  double xiMin    = max( xiMinIn, mXmin * mXmin / s);
  double xiMax    = min( xiMaxIn, sqrt( -tMinIn / SPROTON) );
  if (xiMin > xiMax) return 0.;

  // Integration in xi: check size of affected range and adjust dxi.
  if (xiMax > xiDiv) {
    double xiMinP = max( xiDiv, xiMin);
    int    nxi    = 2 + (xiMax - xiMinP) / dxiRaw;
    double dxi    = (xiMax - xiMinP) / nxi;
    for (int ixi = 0; ixi < nxi; ++ixi)
      sig        += dxi * dsigSDdxi( s, xiMinP + dxi * (ixi + 0.5),
                    tMinIn, tMaxIn, iPart);
  }

  // Integration in ln(xi): check size of affected range and adjust dlnxi.
  if (xiMin < xiDiv) {
    double xiMaxP = min( xiDiv, xiMax);
    int    nlnxi  = 2 + log( xiMaxP / xiMin) / dlnxiRaw;
    double dlnxi  = log( xiMaxP / xiMin) / nlnxi;
    for (int ilnxi = 0; ilnxi < nlnxi; ++ilnxi) {
      double xi   = xiMin * exp( dlnxi * (ilnxi + 0.5));
      sig        += dlnxi * xi * dsigSDdxi( s, xi, tMinIn, tMaxIn, iPart);
    }
  }

  // Done.
  return sig;
}

//==========================================================================

int main() {

  // Energy scale for the fixed-energy plots.
  double eCM   = 10000.;
  double s     = eCM * eCM;

  // Loop through xi for some t values.
  Hist shape1("dsig_SD/dxi, all t",          100, 0., 1.);
  Hist shape2("dsig_SD/dxi, t > -4",         100, 0., 1.);
  Hist shape3("dsig_SD/(dxi dt), t = -0.01", 100, 0., 1.);
  Hist shape4("dsig_SD/(dxi dt), t = -0.1",  100, 0., 1.);
  Hist shape5("dsig_SD/(dxi dt), t = -0.2",  100, 0., 1.);
  Hist shape6("dsig_SD/(dxi dt), t = -0.5",  100, 0., 1.);
  Hist shape7("dsig_SD/(dxi dt), t = -1.0",  100, 0., 1.);
  for (int ixi = 0; ixi < 100; ++ixi) {
    double xi = 0.01 * (ixi + 0.5);
    shape1.fill( xi, dsigSDdxi( s, xi, -s,  0., 0) );
    shape2.fill( xi, dsigSDdxi( s, xi, -4., 0., 0) );
    shape3.fill( xi, dsigSDdtdxi( s, -0.01, xi, 0) );
    shape4.fill( xi, dsigSDdtdxi( s, -0.1,  xi, 0) );
    shape5.fill( xi, dsigSDdtdxi( s, -0.2,  xi, 0) );
    shape6.fill( xi, dsigSDdtdxi( s, -0.5,  xi, 0) );
    shape7.fill( xi, dsigSDdtdxi( s, -1.0,  xi, 0) );
  }
  cout << shape1 << shape2 << shape3 << shape4 << shape5 << shape6 << shape7;

  // Loop through t for some xi values. Note different scale for large xi.
  Hist shape11("dsig_SD/dt, all xi < 0.05",    100, 0., 5.);
  Hist shape12("dsig_SD/(dxi dt), xi = 0.001", 100, 0., 5.);
  Hist shape13("dsig_SD/(dxi dt), xi = 0.01",  100, 0., 5.);
  Hist shape14("dsig_SD/(dxi dt), xi = 0.1",   100, 0., 5.);
  Hist shape15("dsig_SD/(dxi dt), xi = 0.9 (note t scale!)",  100, 0.,  50.);
  Hist shape16("dsig_SD/(dxi dt), xi = 0.99 (note t scale!)", 100, 0., 500.);
  for (int it = 0; it < 100; ++it) {
    double tAbs = 0.05 * (it + 0.5);
    shape11.fill( tAbs, dsigSDdt( s, -tAbs, 0., 0.05, 0) );
    shape12.fill( tAbs, dsigSDdtdxi( s, -tAbs, 0.001, 0) );
    shape13.fill( tAbs, dsigSDdtdxi( s, -tAbs,  0.01, 0) );
    shape14.fill( tAbs, dsigSDdtdxi( s, -tAbs,   0.1, 0) );
    shape15.fill(  10. * tAbs, dsigSDdtdxi( s, -10. * tAbs,    0.9, 0) );
    shape16.fill( 100. * tAbs, dsigSDdtdxi( s, -100. * tAbs,  0.99, 0) );
  }
  cout << shape11 << shape12 << shape13 << shape14 << shape15 << shape16;

  // Loop through xi, in logarithmic xi scale, and split by component.
  Hist shape20("xi * dsig_SD/dxi, lg(xi) scale, all", 100, -10., 0.);
  Hist shape21("xi * dsig_SD/dxi, lg(xi) scale, PPP", 100, -10., 0.);
  Hist shape22("xi * dsig_SD/dxi, lg(xi) scale, PPR", 100, -10., 0.);
  Hist shape23("xi * dsig_SD/dxi, lg(xi) scale, RRP", 100, -10., 0.);
  Hist shape24("xi * dsig_SD/dxi, lg(xi) scale, RRR", 100, -10., 0.);
  Hist shape25("xi * dsig_SD/dxi, lg(xi) scale, pi0", 100, -10., 0.);
  Hist shape26("xi * dsig_SD/dxi, lg(xi) scale, background", 100, -10., 0.);
  Hist shape27("xi * dsig_SD/dxi, lg(xi) scale, resonances", 100, -10., 0.);
  for (int ilgxi = 0; ilgxi < 100; ++ilgxi) {
    double lgxi = -10. + 0.1 * (ilgxi + 0.5);
    double xi   = pow (10., lgxi);
    shape20.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 0) );
    shape21.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 1) );
    shape22.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 2) );
    shape23.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 3) );
    shape24.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 4) );
    shape25.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 5) );
    shape26.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 6) );
    shape27.fill( lgxi, xi * dsigSDdxi( s, xi, -s, 0., 7) );
  }
  cout << shape20 << shape21 << shape22 << shape23 << shape24
       << shape25 << shape26 << shape27;

  // Energy dependence of diffractive cross section.
  Hist shape31("2 * sig_SD(log10(eCM)), 0.01 < xi < 0.05",  80, 1., 5.);
  Hist shape32("2 * sig_SD(log10(eCM)), xi < 0.05, t > -4", 80, 1., 5.);
  Hist shape33("2 * sig_SD(log10(eCM)), xi < 0.05",         80, 1., 5.);
  Hist shape34("2 * sig_SD(log10(eCM)), xi < 0.1",          80, 1., 5.);
  Hist shape35("2 * sig_SD(log10(eCM)), xi < 0.2",          80, 1., 5.);
  Hist shape36("2 * sig_SD(log10(eCM)), xi < 0.4",          80, 1., 5.);
  Hist shape41("2 * sig_SD(log10(eCM)), xi < 0.05, PPP",    80, 1., 5.);
  Hist shape42("2 * sig_SD(log10(eCM)), xi < 0.05, PPR",    80, 1., 5.);
  Hist shape43("2 * sig_SD(log10(eCM)), xi < 0.05, RRP",    80, 1., 5.);
  Hist shape44("2 * sig_SD(log10(eCM)), xi < 0.05, RRR",    80, 1., 5.);
  Hist shape45("2 * sig_SD(log10(eCM)), xi < 0.05, pi0",    80, 1., 5.);
  Hist shape46("2 * sig_SD(log10(eCM)), xi < 0.05, background", 80, 1., 5.);
  Hist shape47("2 * sig_SD(log10(eCM)), xi < 0.05, resonances", 80, 1., 5.);
  for (int iE = 20; iE < 100; ++iE) {
    double eLog   = 0.05 * (iE + 0.5);
    double eCMnow = pow( 10., eLog);
    double sNow   = eCMnow * eCMnow;
    shape31.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.01, 0.05, 0) );
    shape32.fill ( eLog, 2. * sigSD( sNow, -4.,0., 0.0,  0.05, 0) );
    shape33.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 0) );
    shape34.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.10, 0) );
    shape35.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.20, 0) );
    shape36.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.40, 0) );
    shape41.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 1) );
    shape42.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 2) );
    shape43.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 3) );
    shape44.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 4) );
    shape45.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 5) );
    shape46.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 6) );
    shape47.fill ( eLog, 2. * sigSD( sNow, -s, 0., 0.0,  0.05, 7) );
  }
  cout << shape31 << shape32 << shape33 << shape34 << shape35 << shape36
       << shape41 << shape42 << shape43 << shape44 << shape45
       << shape46 << shape47;

  // Low-mass resonance contributions.
  Hist shape51("dsig/dMX mX < 20", 100, 0., 20.);
  Hist shape52("dsig/dMX mX < 4, resonances only", 100, 0., 4.);
  for (int iMX = 0; iMX < 100; ++iMX) {
    double mXnow = 0.2 * (iMX + 0.5);
    double xiNow = mXnow * mXnow / s;
    shape51.fill( mXnow, dsigSDdxi( s, xiNow) * 2. * mXnow / s);
    mXnow *= 0.2;
    xiNow = mXnow * mXnow / s;
    shape52.fill( mXnow, dsigSDdxi( s, xiNow, -s, 0., 7) * 2. * mXnow / s);
  }
  cout << shape51 << shape52;

  // Figure 8 in ABMST.
  Hist fig08("Fig 8: dsig/dtdxi, eCM = 17.58, t = -0.131" , 100, 0., 0.2);
  double s08 = pow2(17.58);
  double t08 = -0.131;
  for (int ixi = 0; ixi < 100; ++ixi) {
    double xi = 0.002 * (ixi + 0.5);
    double dsig = dsigSDdtdxi( s08, t08, xi, 0);
    fig08.fill( xi, dsig);
  }
  cout << fig08;

  // Figure 9 in ABMST.
  Hist fig09("Fig 9: dsig/dtdxi, eCM = 53.67, t = -0.52" , 100, 0., 0.2);
  double s09 = pow2(53.67);
  double t09 = -0.52;
  for (int ixi = 0; ixi < 100; ++ixi) {
    double xi = 0.002 * (ixi + 0.5);
    double dsig = dsigSDdtdxi( s09, t09, xi, 0);
    fig09.fill( xi, dsig);
  }
  cout << fig09;

  // Figure 10 in ABMST.
  Hist fig10("Fig 10: dsig/dMX2dt, eCM = 23.7, t = -0.05" , 100, 0., 10.);
  double s10 = pow2(23.7);
  double t10 = -0.05;
  for (int iMX2 = 0; iMX2 < 100; ++iMX2) {
    double mX2 = 0.1 * (iMX2 + 0.5);
    double xi  = mX2 / s10;
    double dsig = dsigSDdtdxi( s10, t10, xi, 0);
    fig10.fill( mX2, dsig / s10);
  }
  cout << fig10;

  // Figure 11 in ABMST.
  Hist fig11("Fig 11: dsig/dMX2dt(log10(MX2)), eCM = 23.7, t = -0.05" ,
    100, 0., 2.);
  double s11 = pow2(23.7);
  double t11 = -0.05;
  for (int iMX2 = 0; iMX2 < 100; ++iMX2) {
    double lg10mX2 = 0.02 * (iMX2 + 0.5);
    double mX2 = pow(10., lg10mX2);
    double xi  = mX2 / s10;
    double dsig = dsigSDdtdxi( s11, t11, xi, 0);
    fig11.fill( lg10mX2, dsig / s11);
  }
  cout << fig11;

  // Figure 12 in ABMST.
  Hist fig12("Fig 12: dsig/dt, eCM = 30.5, xi < 0.05 (?)" , 100, 0., 2.);
  double s12 = pow2(30.5);
  for (int it = -0; it < 100; ++it) {
    double tAbs = 0.02 * (it + 0.5);
    double dsig = dsigSDdt( s12, -tAbs, 0.0, 0.05, 0);
    fig12.fill( tAbs, dsig);
  }
  cout << fig12;

  // Figure 13 in ABMST.
  Hist fig13("Fig 13: dsig/dt, eCM = 38.3, xi < 0.05 (?)" , 100, 0., 2.);
  double s13 = pow2(38.3);
  for (int it = -0; it < 100; ++it) {
    double tAbs = 0.02 * (it + 0.5);
    double dsig = dsigSDdt( s13, -tAbs, 0.0, 0.05, 0);
    fig13.fill( tAbs, dsig);
  }
  cout << fig13;

  // Figure 14 in ABMST.
  Hist fig14("Fig 14: dsig/dt, eCM = 546, xi < 0.05 (?)" , 100, 0., 2.);
  double s14 = pow2(546.);
  for (int it = -0; it < 100; ++it) {
    double tAbs = 0.02 * (it + 0.5);
    double dsig = dsigSDdt( s14, -tAbs, 0.0, 0.05, 0);
    fig14.fill( tAbs, dsig);
  }
  cout << fig14;

  // Figure 15 in ABMST: table of integrated SD cross section.
  double ecmTab[8] = { 14., 20., 30., 40., 60., 100., 200., 500.};
  double sigABMST[8] = { 3.5, 4.5, 5.2, 5.6, 6.5, 6.8, 7.8, 9.3};
  cout << "\n     eCM    code   ABMST " << endl;
  for (int iTab = 0; iTab < 8; ++iTab) {
    double sTab   = pow2(ecmTab[iTab]);
    double sigTab = 2. * sigSD( sTab, -sTab, 0.0, 0.0,  0.05, 0);
    cout << fixed << setprecision(0) << setw(8) << ecmTab[iTab]
         << setprecision(2) << setw(8) << sigTab << setw(8)
         << sigABMST[iTab] << endl;
  }

  return 0;
}
