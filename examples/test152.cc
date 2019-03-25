// main51.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of LHAPDF interface and whether PDF's behave sensibly.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

class LHAGrid1 : public PDF {

public:

  // Constructor.
  LHAGrid1(int idBeamIn = 2212, int iFitIn = 1,
    string xmlPath = "../share/Pythia8/xmldoc/", Info* infoPtr = 0)
    : PDF(idBeamIn), doExtraPol(false), pdfGrid(NULL), pdfSlope(NULL) {
    init( iFitIn, xmlPath, infoPtr); };

  // Destructor.
  ~LHAGrid1() { if (pdfGrid) { for (int iid = 0; iid < 12; ++iid) {
    for (int ix = 0; ix < nx; ++ix) delete[] pdfGrid[iid][ix];
    delete[] pdfGrid[iid]; } delete[] pdfGrid; }
    if (pdfSlope) { for (int iid = 0; iid < 12; ++iid) delete[] pdfSlope[iid];
    delete[] pdfSlope;} };

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool doExtraPolIn) {doExtraPol = doExtraPolIn;}

private:

  // Variables to be set during code initialization.
  bool   doExtraPol;
  int    iFit, nx, nq;
  double xMin, xMax, qMin, qMax, pdfVal[12];
  vector<double> xGrid, lnxGrid, qGrid, lnqGrid;
  double*** pdfGrid;
  double** pdfSlope;

  // Initialization of data array.
  void init( int iFitIn, string xmlPath, Info* infoPtr);

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // Interpolation in the grid for a given PDF flavour.
  void xfxevolve(double x, double Q2);

};


//--------------------------------------------------------------------------

// Initialize PDF: read in data grid from file.

void LHAGrid1::init(int iFitIn, string xmlPath, Info* infoPtr) {

  // Choice of fit among possibilities. Some local variables.
  iFit = iFitIn;
  string line;
  double xNow, qNow, pdfNow;
  vector<int> idGridMap;
  int idNow, idNowMap;

  // Open files from which grids should be read in.
  if (xmlPath[ xmlPath.length() - 1 ] != '/') xmlPath += "/";
  string         dataFile = "NNPDF23_lo_as_0119_qed_0000.dat";
  //if (iFit == 1) dataFile = "pomH1FitA.data";
  //if (iFit == 2) dataFile = "pomH1FitB.data";
  ifstream is( (xmlPath + dataFile).c_str() );
  if (!is.good()) {
    if (infoPtr != 0) infoPtr->errorMsg("Error from LHAGrid1::init: "
      "the requested parametrization file was not found");
    else cout << " Error from LHAGrid1::init: "
      << "the requested parametrization file was not found" << endl;
    isSet = false;
    return;
  }

  // Skip ahead and then read in x and Q grids.
  for (int i = 0; i < 3; ++i) getline( is, line);
  getline( is, line);
  istringstream isx(line);
  while (isx >> xNow) {
    xGrid.push_back( xNow);
    lnxGrid.push_back( log(xNow));
  }
  nx   = xGrid.size();
  xMin = xGrid.front();
  xMax = xGrid.back();
  getline( is, line);
  istringstream isq(line);
  while (isq >> qNow) {
    qGrid.push_back( qNow);
    lnqGrid.push_back( log(qNow));
  }
  nq   = qGrid.size();
  qMin = qGrid.front();
  qMax = qGrid.back();
  //for (int ix = 0; ix < nx; ++ix) cout << ix << " " << scientific
  //  << xGrid[ix] << endl;
  //for (int iq = 0; iq < nq; ++iq) cout << iq << " " << scientific
  //  << qGrid[iq] << endl;

  // Read in flavour grid and decide flavour mapping.
  getline( is, line);
  istringstream isid(line);
  while (isid >> idNow) {
    idNowMap = -1;
    if (idNow == 21 || idNow == 0) idNowMap = 0;
    if (idNow > 0 && idNow < 6) idNowMap = idNow;
    if (idNow < 0 && idNow > -6) idNowMap = 5 - idNow;
    if (idNow == 22) idNowMap = 11;
    idGridMap.push_back( idNowMap);
    //cout << setw(4) << idNow << setw(5) << idNowMap << endl;
  }
  int nid = idGridMap.size();

  // Create array big enough to hold (flavour, x, Q) grid.
  pdfGrid = new double**[12];
  for (int iid = 0; iid < 12; ++iid) {
    pdfGrid[iid] = new double*[nx];
    for (int ix = 0; ix < nx; ++ix) {
      pdfGrid[iid][ix] = new double[nq];
      for (int iq = 0; iq < nq; ++iq) pdfGrid[iid][ix][iq] = 0.;
    }
  }

  // Read in data grid, line by line.
  for (int ix = 0; ix < nx; ++ix)
  for (int iq = 0; iq < nq; ++iq) {
    getline( is, line);
    istringstream ispdf(line);
    for (int iid = 0; iid < nid; ++iid) {
      ispdf >> pdfNow;
      if (idGridMap[iid] >= 0) pdfGrid[idGridMap[iid]][ix][iq] = pdfNow;
    }
  }

  // For extrapolation to small x: create array for b values of x^b shape.
  pdfSlope = new double*[12];
  for (int iid = 0; iid < 12; ++iid) {
    pdfSlope[iid] = new double[nq];
    for (int iq = 0; iq < nq; ++iq) { pdfSlope[iid][iq] =
      ( min( pdfGrid[iid][0][iq], pdfGrid[iid][1][iq]) > 1e-5)
      ? ( log(pdfGrid[iid][1][iq]) - log(pdfGrid[iid][0][iq]) )
      / (lnxGrid[1] - lnxGrid[0]) : 0.;
    }
  }

}

//--------------------------------------------------------------------------

void LHAGrid1::xfUpdate(int , double x, double Q2) {

  // Update using NNPDF routine, within allowed (x, q) range.
  xfxevolve(x,Q2);

  // Then transfer to Pythia8 notation.
  xg     = pdfVal[0];
  xu     = pdfVal[2];
  xd     = pdfVal[1];
  xubar  = pdfVal[7];
  xdbar  = pdfVal[6];
  xs     = pdfVal[3];
  xsbar  = pdfVal[8];
  xc     = 0.5 * (pdfVal[4] + pdfVal[9]);
  xb     = 0.5 * (pdfVal[5] + pdfVal[10]);
  xgamma = pdfVal[11];

  // Subdivision of valence and sea.
  xuVal  = xu - xubar;
  xuSea  = xubar;
  xdVal  = xd - xdbar;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav  = 9;

}

//--------------------------------------------------------------------------

void LHAGrid1::xfxevolve(double x, double Q2) {

  // Find if (x, Q) inside our outside grid.
  double q = sqrt(Q2);
  int inx  = (x <= xMin) ? -1 : ((x >= xMax) ? 1 : 0);
  int inq  = (q <= qMin) ? -1 : ((q >= qMax) ? 1 : 0);

  // Find grid values surrounding x; freeze at border.
  int    minx = 0;
  int    maxx = nx -1;
  double dlnx = 0.5;
  if (inx == 0) {
    int midx;
    while (maxx - minx > 1) {
      midx = (minx + maxx) / 2;
      if (x < xGrid[midx]) maxx = midx;
      else                 minx = midx;
    }
    dlnx = (log(x) - lnxGrid[minx]) / ( lnxGrid[maxx] - lnxGrid[minx]);
  } else if (inx == -1) maxx = minx;
    else if (inx ==  1) minx = maxx;

  // Find grid values surrounding q; freeze at border.
  int    minq = 0;
  int    maxq = nq - 1;
  double dlnq = 0.5;
  if (inq == 0) {
    int midq;
    while (maxq - minq > 1) {
      midq = (minq + maxq) / 2;
      if (q < qGrid[midq]) maxq = midq;
      else                 minq = midq;
    }
    dlnq = (log(q) - lnqGrid[minq]) / ( lnqGrid[maxq] - lnqGrid[minq]);
  } else if (inq == -1) maxq = minq;
    else if (inq ==  1) minq = maxq;

  // Special: extrapolate to small x.
  if (inx == -1 && doExtraPol) {
    for (int iid = 0; iid < 12; ++iid)
      pdfVal[iid] = (1. - dlnq) * pdfGrid[iid][0][minq]
                  * pow( x / xMin, pdfSlope[iid][minq])
                  +       dlnq  * pdfGrid[iid][0][maxq]
                  * pow( x / xMin, pdfSlope[iid][maxq]);

  // Normal: interpolate between grid elements.
  } else {
    for (int iid = 0; iid < 12; ++iid)
      pdfVal[iid] = (1. - dlnx) * (1. - dlnq) * pdfGrid[iid][minx][minq]
                  +       dlnx  * (1. - dlnq) * pdfGrid[iid][maxx][minq]
                  + (1. - dlnx) *       dlnq  * pdfGrid[iid][minx][maxq]
                  +       dlnx  *       dlnq  * pdfGrid[iid][maxx][maxq];
  }

}

//==========================================================================

// Integration to check momentum sum rule.

double integrate(PDF* nowPDF, double Q2) {

  // Number of points, x ranges and initial values.
  int    nLin  = 980;
  int    nLog  = 1000;
  double xLin  = 0.02;
  double xLog  = 1e-8;
  double dxLin = (1. - xLin) / nLin;
  double dxLog = log(xLin / xLog) / nLog;
  double sum   = 0.;
  double x, sumNow;

  // Integration at large x in linear steps.
  for (int iLin = 0; iLin < nLin; ++iLin) {
    x      = xLin + (iLin + 0.5) * dxLin;
    sumNow = nowPDF->xf( 21, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLin * sumNow;
  }

  // Integration at small x in logarithmic steps.
  for (int iLog = 0; iLog < nLog; ++iLog) {
    x      = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
    sumNow = nowPDF->xf( 21, x, Q2);
    for (int i = 1; i < 6; ++i)
      sumNow += nowPDF->xf( i, x, Q2) + nowPDF->xf( -i, x, Q2);
    sum   += dxLog * x * sumNow;
  }

  // Done.
  return sum;

}

//==========================================================================

int main() {

  // The Pythia class itself is not used, but some facilities that come along.
  //Pythia pythia;

  // Chosen new PDF set; LHAPDF5 file name conventions.
  //string pdfSet = "LHAPDF5:cteq5l.LHgrid";
  //string pdfSet = "LHAPDF5:cteq61.LHpdf";
  //string pdfSet = "LHAPDF5:cteq61.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2004nlo.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2001lo.LHgrid";
  //string pdfSet = "LHAPDF5:MRST2007lomod.LHgrid";

  // Pointers to old default and new tryout PDF sets.
  //Info info;
  PDF* oldPDF = new NNPDF(2212, 2);
  PDF* newPDF = new LHAGrid1(2212, 1);

  // Alternative: compare two Pomeron PDF's. Boost second by factor 2.
  //PDF* oldPDF = new PomFix( 990, -0.2, 2.5, 0., 3., 0.4, 0.5);
  //PDF* newPDF = new PomH1Jets( 990, 2.);
  //PDF* oldPDF = new PomH1FitAB( 990, 2);
  //PDF* newPDF = new PomH1FitAB( 990, 3);

  // Allow extrapolation of PDF's beyond x and Q2 boundaries, at own risk.
  // Default behaviour is to freeze PDF's at boundaries.
  //newPDF->setExtrapolate(true);

  // Histogram F(x, Q2) = (9/4) x*g(x, Q2) + sum_{i = q, qbar} x*f_i(x, Q2)
  // for range 10^{-8} < x < 1 logarithmic in x and for Q2 = 4 and 100.
  Hist oldF4("F( x, Q2 = 4) old", 80 , -8., 0.);
  Hist newF4("F( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratF4("F( x, Q2 = 4) new/old", 80 , -8., 0.);
  Hist oldF100("F( x, Q2 = 100) old", 80 , -8., 0.);
  Hist newF100("F( x, Q2 = 100) new", 80 , -8., 0.);
  Hist ratF100("F( x, Q2 = 100) new/old", 80 , -8., 0.);

  // Loop over the two Q2 values.
  for (int iQ = 0; iQ < 2; ++iQ) {
    double Q2 = (iQ == 0) ? 4. : 100;

    // Loop over x values, in a logarithmic scale
    for (int iX = 0; iX < 80; ++iX) {
      double xLog = -(0.1 * iX + 0.05);
      double x = pow( 10., xLog);

      // Evaluate old summed PDF.
      double oldSum = 2.25 * oldPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        oldSum += oldPDF->xf( i, x, Q2) + oldPDF->xf( -i, x, Q2);
      if (iQ == 0) oldF4.fill ( xLog, oldSum );
      else       oldF100.fill ( xLog, oldSum );

      // Evaluate new summed PDF.
      double newSum = 2.25 * newPDF->xf( 21, x, Q2);
      for (int i = 1; i < 6; ++i)
        newSum += newPDF->xf( i, x, Q2) + newPDF->xf( -i, x, Q2);
      if (iQ == 0) newF4.fill ( xLog, newSum );
      else       newF100.fill ( xLog, newSum );

    //End loops over x and Q2 values.
    }
  }

  // Show F(x, Q2) and their ratio new/old.
  ratF4 = newF4 / oldF4;
  ratF100 = newF100 / oldF100;
  cout << oldF4 << newF4 << ratF4 << oldF100 << newF100 << ratF100;

  // Histogram momentum sum as a function of Q2 (or rather log10(Q2)).
  Hist oldXSum("momentum sum(log10(Q2)) old", 100, -2., 8.);
  Hist newXSum("momentum sum(log10(Q2)) new", 100, -2., 8.);

  // Loop over Q2 values.
  for (int iQ = 0; iQ < 100; ++iQ) {
    double log10Q2 = -2.0 + 0.1 * iQ + 0.05;
    double Q2 = pow( 10., log10Q2);

    // Evaluate old and new momentum sums.
    double oldSum = integrate( oldPDF, Q2);
    oldXSum.fill( log10Q2, oldSum);
    double newSum = integrate( newPDF, Q2);
    newXSum.fill( log10Q2, newSum);
  }

  // Show momentum sum as a function of Q2.
  cout << oldXSum << newXSum;

  // Pythia instance to get random numbers.
  Pythia pythia;

  // Comparison of value at random points.
  Hist diff("(new - old)/(new + old)", 100, -0.02, 0.02);
  int id[6] = { 21, 1, -1, 2, -2, 3};
  double x, Q2, oldVal, newVal, diffrat;
  for (int iTry = 0; iTry < 10000000; ++iTry) {
    x  = pow( 1e-6, pythia.rndm.flat() );
    Q2 = pow( 1e6,  pythia.rndm.flat() );
    for (int iid = 0 ; iid < 6; ++iid) {
      oldVal  = oldPDF->xf( id[iid], x, Q2);
      newVal  = newPDF->xf( id[iid], x, Q2);
      diffrat = (newVal - oldVal) / (newVal + oldVal);
      if (oldVal > 1e-3 || newVal > 1e-3) {
        diff.fill( diffrat );
        if (abs(diffrat) > 0.05) cout << scientific << " x = " << x << " Q2 = "
          << Q2 << " id = " << id[iid] << " old = " << oldVal << " new = "
          << newVal << " ratio = " <<diffrat << endl;
      }
    }
  }
  cout << diff;

  // Done.
  delete oldPDF;
  delete newPDF;
  return 0;
}
