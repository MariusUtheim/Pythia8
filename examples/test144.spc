# SOFTSUSY SUGRA calculation
# SOFTSUSY3.6.2 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.6.2       # version number
Block MODSEL  # Select model
     1    1   # sugra
Block SMINPUTS             # Standard Model inputs
     1    1.27916000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.18400000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.18000000e+00   # mb(mb)
     6    1.73500000e+02   # Mtop(pole)
     7    1.77699000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    3.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    2.88700000e+03   # m0
     2    5.58000000e+02   # m12
     5   -5.77400000e+03   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.76304997e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=1 Desired accuracy=1.00000000e-04 Achieved accuracy=2.21344586e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03881004e+01   # MW
        25     1.25919331e+02   # h0
        35     2.80496467e+03   # H0
        36     2.80358839e+03   # A0
        37     2.80536824e+03   # H+
   1000021     1.39974387e+03   # ~g
   1000022     2.46356232e+02   # ~neutralino(1)
   1000023     4.80304280e+02   # ~neutralino(2)
   1000024     4.80465221e+02   # ~chargino(1)
   1000025    -2.18220567e+03   # ~neutralino(3)
   1000035     2.18343928e+03   # ~neutralino(4)
   1000037     2.18414506e+03   # ~chargino(2)
   1000001     3.06134477e+03   # ~d_L
   1000002     3.06045108e+03   # ~u_L
   1000003     3.06119490e+03   # ~s_L
   1000004     3.06030116e+03   # ~c_L
   1000005     1.96517713e+03   # ~b_1
   1000006     8.55137645e+02   # ~t_1
   1000011     2.90312431e+03   # ~e_L
   1000012     2.90172435e+03   # ~nue_L
   1000013     2.90227501e+03   # ~mu_L
   1000014     2.90087473e+03   # ~numu_L
   1000015     2.35215740e+03   # ~stau_1
   1000016     2.64832721e+03   # ~nu_tau_L
   2000001     3.05003168e+03   # ~d_R
   2000002     3.05171202e+03   # ~u_R
   2000003     3.04975840e+03   # ~s_R
   2000004     3.05168282e+03   # ~c_R
   2000005     2.61502079e+03   # ~b_2
   2000006     1.98733103e+03   # ~t_2
   2000011     2.89053666e+03   # ~e_R
   2000013     2.88882461e+03   # ~mu_R
   2000015     2.65179388e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -3.43366896e-02       # alpha - evaluated at p^2=0
Block nmix                  # neutralino mixing matrix
  1  1     9.99765148e-01   # N_{1,1}
  1  2    -1.68545796e-03   # N_{1,2}
  1  3     2.13372254e-02   # N_{1,3}
  1  4    -3.39566271e-03   # N_{1,4}
  2  1     2.52273729e-03   # N_{2,1}
  2  2     9.99237437e-01   # N_{2,2}
  2  3    -3.77162105e-02   # N_{2,3}
  2  4     9.78095762e-03   # N_{2,4}
  3  1    -1.26447999e-02   # N_{3,1}
  3  2     1.97795046e-02   # N_{3,2}
  3  3     7.06590996e-01   # N_{3,3}
  3  4     7.07232667e-01   # N_{3,4}
  4  1    -1.74181892e-02   # N_{4,1}
  4  2     3.36225267e-02   # N_{4,2}
  4  3     7.06294113e-01   # N_{4,3}
  4  4    -7.06905056e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.98574258e-01   # U_{1,1}
  1  2    -5.33802600e-02   # U_{1,2}
  2  1     5.33802600e-02   # U_{2,1}
  2  2     9.98574258e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.99902983e-01   # V_{1,1}
  1  2    -1.39292962e-02   # V_{1,2}
  2  1     1.39292962e-02   # V_{2,1}
  2  2     9.99902983e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     1.40928212e-01   # F_{11}
  1  2     9.90019818e-01   # F_{12}
  2  1     9.90019818e-01   # F_{21}
  2  2    -1.40928212e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.98737680e-01   # F_{11}
  1  2     5.02299461e-02   # F_{12}
  2  1    -5.02299461e-02   # F_{21}
  2  2     9.98737680e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     8.37273577e-02   # F_{11}
  1  2     9.96488700e-01   # F_{12}
  2  1     9.96488700e-01   # F_{21}
  2  2    -8.37273577e-02   # F_{22}
Block gauge Q= 1.23371981e+03  # SM gauge couplings
     1     3.61851076e-01   # g'(Q)MSSM DRbar
     2     6.37501377e-01   # g(Q)MSSM DRbar
     3     1.04356134e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.23371981e+03  
  3  3     8.40270174e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.23371981e+03  
  3  3     3.79060805e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.23371981e+03  
  3  3     3.03452513e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.23371981e+03 # Higgs mixing parameters
     1     2.18541188e+03    # mu(Q)MSSM DRbar
     2     2.91704051e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43371952e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.03327733e+07    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.23371981e+03  # MSSM DRbar SUSY breaking parameters
     1     2.44873449e+02      # M_1(Q)
     2     4.50547622e+02      # M_2(Q)
     3     1.23344411e+03      # M_3(Q)
    21     3.07901797e+06      # mH1^2(Q)
    22    -4.79627971e+06      # mH2^2(Q)
    31     2.89703836e+03      # meL(Q)
    32     2.89619051e+03      # mmuL(Q)
    33     2.63993664e+03      # mtauL(Q)
    34     2.88887447e+03      # meR(Q)
    35     2.88716416e+03      # mmuR(Q)
    36     2.34053340e+03      # mtauR(Q)
    41     3.03328803e+03      # mqL1(Q)
    42     3.03313896e+03      # mqL2(Q)
    43     1.93676552e+03      # mqL3(Q)
    44     3.02740713e+03      # muR(Q)
    45     3.02737803e+03      # mcR(Q)
    46     8.03759692e+02      # mtR(Q)
    47     3.02707513e+03      # mdR(Q)
    48     3.02680288e+03      # msR(Q)
    49     2.57782567e+03      # mbR(Q)
Block au Q= 1.23371981e+03  
  1  1    -5.17118642e+03      # Au(Q)MSSM DRbar
  2  2    -5.17107244e+03      # Ac(Q)MSSM DRbar
  3  3    -2.95897895e+03      # At(Q)MSSM DRbar
Block ad Q= 1.23371981e+03  
  1  1    -6.69300889e+03      # Ad(Q)MSSM DRbar
  2  2    -6.69274943e+03      # As(Q)MSSM DRbar
  3  3    -5.57878035e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.23371981e+03  
  1  1    -5.51413981e+03      # Ae(Q)MSSM DRbar
  2  2    -5.51260819e+03      # Amu(Q)MSSM DRbar
  3  3    -5.04693729e+03      # Atau(Q)MSSM DRbar
