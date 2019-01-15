//#include <iostream>
//#include "Pythia8/ResonanceData.h"
//#include "../tests.h"
//
//using namespace Pythia8;
//
//void test_ClebschGordan() {
//  cout << "Testing Clebsch-Gordan relations..." << endl;
//
//  ResonanceData r;
////printf("%.3f\n", r.getClebschGordan2(2, -2, 3, 1, 1, -1));
////return;
//
//  // Test JM completeness
//  for (int j1 = 0; j1 <= 3; ++j1)
//  for (int j2 = 0; j2 <= 3 && j1 + j2 < 6; ++j2)
//  for (int m1 = -j1; m1 <= j1; m1 += 2)
//  for (int m2 = -j2; m2 <= j2; m2 += 2) {
//    double sum = 0;
//    for (int j = abs(j1 - j2); j <= j1 + j2; j += 2)
//    for (int m = -j; m <= j; m += 2)
//      if (m == m1 + m2)
//        sum += r.getClebschGordan2(j1, m1, j2, m2, j, m);
//    if (abs(sum - 1.) > 10e-12)
//      printf(" ERROR: sum of JM states is not 1 for j1=%d m1=%d j2=%d m2=%d (sum=%.16f)\n", j1, m1, j2, m2, sum);
//  }
//
//  // Test m1m2 completeness
//  for (int j1 = 0; j1 <= 3; ++j1)
//  for (int j2 = 0; j2 <= 3 && j1 + j2 < 6; ++j2)
//  for (int j = abs(j1 - j2); j <= j1 + j2; j += 2)
//  for (int m = -j; m <= j; m += 2) {
//    double sum = 0;
//    for (int m1 = -j1; m1 <= j1; m1 += 2)
//    for (int m2 = -j2; m2 <= j2; m2 += 2)
//      if (m1 + m2 == m)
//        sum += r.getClebschGordan2(j1, m1, j2, m2, j, m);
//    if (abs(sum - 1.) > 10e-12)
//      printf(" ERROR: sum over m1,m2 is not 1 for j1=%d j2=%d J=%d M=%d (sum=%.16f)\n", j1, j2, j, m, sum);
//  }
//
//  // Test symmetries
//  for (int j1 = 0; j1 <= 3; ++j1)
//  for (int m1 = -j1; m1 <= j1; m1 += 2)
//  for (int j2 = 0; j2 <= 3 && j1 + j2 < 6; ++j2)
//  for (int m2 = -j2; m2 <= j2; m2 += 2)
//  for (int j = abs(j1 - j2); j <= j1 + j2; j += 2)
//  for (int m = -j; m <= j; m += 2) 
//  if (m1 + m2 == m) {
//    if (r.getClebschGordan2(j1, -m1, j2, -m2, j, -m) != r.getClebschGordan2(j1, m1, j2, m2, j, m))
//      cout << " ERROR: CG(j1,-m1,j2,-m2,j,-m) != CG(j1,m1,j2,m2,j,m)" << endl;
//    if (r.getClebschGordan2(j2, m2, j1, m1, j, m) != r.getClebschGordan2(j1, m1, j2, m2, j, m))
//      cout << " ERROR: CG(j2,m2,j1,m,j,-m) != CG(j1,m1,j2,m2,j,m)" << endl;
//
//    if (j <= 3 && j + j1 < 6) {
//      double lhs, rhs;
//      
//      // Note: the j's here are twice that the j's that are usually used,
//      // e.g. for a spin 1/2 particle, j = 1 here.
//      lhs = r.getClebschGordan2(j1, m1, j, -m, j2, -m2);
//      rhs = (j2 + 1.) / (j + 1.) * r.getClebschGordan2(j1, m1, j2, m2, j, m);
//      if (lhs != rhs)
//        printf(" ERROR: CG(j1,m1,j,-m,j2,-m2) != (2j2+1)/(2j+1) CG(j1,m1,j2,m2,j,m)\n\t\t for j1=%d m1=%d j2=%d m2=%d J=%d M=%d (lhs=%f, rhs=%f, factor=%f)\n\n",
//          j1, m1, j2, m2, j, m, lhs, rhs, (j2 + 1.) / (j + 1.));
//
//      lhs = r.getClebschGordan2(j, m, j1, -m1, j2, m2);
//      rhs = (j2 + 1.) / (j + 1.) * r.getClebschGordan2(j1, m1, j2, m2, j, m);
//      if (lhs != rhs)
//        printf(" ERROR: CG(j,m,j1,-m1,j2,m2) != (2j2+1)/(2j+1) CG(j1,m1,j2,m2,j,m)\n\t\t for j1=%d m1=%d j2=%d m2=%d J=%d M=%d (lhs=%f, rhs=%f, factor=%f)\n\n",
//          j1, m1, j2, m2, j, m, lhs, rhs, (j2 + 1.) / (j + 1.));
//
//    }
//  }
//
//  cout << "test complete." << endl << endl;
//}
//