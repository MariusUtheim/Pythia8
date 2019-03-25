// test250.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Class showing how to extract the PDF of a second hard interaction,
// given the nature of the first interaction. 
// In this example used to check the momentum sum rule.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Evaluate the proton PDF of a second collision, given a first collision.
//
// Select the ordinary PDF set in the constructor of the SecondPDF class;
// see the PYTHIA HTML manual, page "PDF Selection", word "PDF:pSet".
// The xf2 method then returns the xf_i(x, Q^2) value at (id_2, x_2, Q_2^2) 
// given that a first has already been extracted at (id_1, x_1, Q_12).
//
// By comparison, the xf1 method returns corresponding normal PDFs.
// The product of xf1 * xf2 gives the doubly differential distribution.
// Symmetrize if the two interactions 1 and 2 happen at comparable scales:
// ( xf1( 1) * xf2( 2; 1) + xf1( 2) * xf2( 1; 2) ) / 2.
//
// Warning: the evaluation of the second PDF is different whether the 
// first picks up a valence quark or not, and the outcome can therefore 
// be either of two values, mixed in appropriate fractions. An average of 
// several evaluations thus is needed for a smooth PDF curve, but averaging 
// occurs automatically in the context of Monte Carlo event generation.  

class SecondPDF {

public:

  // Constructor, giving PDF selection as input.
  SecondPDF( string pdfSet) {
    // Set up fictitious pp collision process.  
    pythia.readString("HardQCD:hardbbbar = on");
    pythia.readString("PartonLevel:MPI = off");
    // Select appropriate PDF set and initialize. 
    pythia.settings.word("PDF:pSet", pdfSet);
    pythia.init();
  }

  // Normal PDF, as for first interaction. (For simple comparisons.)
  double xf1( int id, double x, double Q2) {
    return pythia.beamA.xf( id, x, Q2);
  }
    
  // PDF for second interaction, given a first one.    
  double xf2( int id2, double x2, double Q22, int id1, double x1, double Q12) {
    // Set one incoming beam as if first parton has been extracted from it.
    pythia.beamA.clear();
    pythia.beamA.append( 0, id1, x1);
    pythia.beamA.xfISR(  0, id1, x1, Q12);
    pythia.beamA.pickValSeaComp();
    // Evaluate PDF of second parton given the first one.
    return pythia.beamA.xfMPI( id2, x2, Q22);
  }

private:

  // Pythia instance.
  Pythia pythia;

};

//==========================================================================

// Main program does simple integration of second PDF to check momentum sum. 

int main() {

  // Create object for PDF evaluation. Argument gives PDF set, see above.
  SecondPDF secPDF("8");

  // Number of integration points, x ranges and local variables.
  int    nLin  = 980;
  int    nLog  = 1000;
  double xLin  = 0.02;
  double xLog  = 1e-8;
  double dxLin = (1. - xLin) / nLin;
  double dxLog = log(xLin / xLog) / nLog;
  double x2, sumNow;

  // Set Q2 scales for the two interactions.
  double Q12 = 100.;
  double Q22 = 20.;

  // Loop over flavour and x values for the first interaction.
  for (int iid1 = -5; iid1 < 6; ++iid1)
  for (int ix1 = 0; ix1 < 7; ++ ix1) {
    int id1    = (iid1 == 0) ? 21 : iid1;
    double x1  = (ix1 < 4) ? pow( 10., ix1 - 4) : 0.1 * (ix1 - 2);

    // Initialize momentum sum to begin integration.
    double sum = 0.;

    // Integration at large x in linear steps.
    for (int iLin = 0; iLin < nLin; ++iLin) {
      x2     = xLin + (iLin + 0.5) * dxLin;
      sumNow = secPDF.xf2( 21, x2, Q22, id1, x1, Q12) 
             + secPDF.xf2( 22, x2, Q22, id1, x1, Q12);
      for (int i = 1; i < 6; ++i)
        sumNow += secPDF.xf2(  i, x2, Q22, id1, x1, Q12) 
                + secPDF.xf2( -i, x2, Q22, id1, x1, Q12);
      sum   += dxLin * sumNow;
    }

    // Integration at small x in logarithmic steps.
    for (int iLog = 0; iLog < nLog; ++iLog) {
      x2     = xLog * pow( xLin / xLog, (iLog + 0.5) / nLog );
      sumNow = secPDF.xf2( 21, x2, Q22, id1, x1, Q12) 
             + secPDF.xf2( 22, x2, Q22, id1, x1, Q12);
      for (int i = 1; i < 6; ++i)
        sumNow += secPDF.xf2(  i, x2, Q22, id1, x1, Q12) 
                + secPDF.xf2( -i, x2, Q22, id1, x1, Q12);
      sum   += dxLog * x2 * sumNow;
    }

    // Result.
    cout << " For id1 = " << setw(2) << id1 << " and x1 = " << fixed
         << setprecision(5) << x1 << " the momentum sum is "
         << setprecision(5) << x1 + sum << endl;
  }

  // Done.
  return 0;
}
