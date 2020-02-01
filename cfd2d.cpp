#include "cfd2d.h"
#include <math.h>
#include <cstdio>
#include <stdlib.h>
//#include <conio.h>

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
////   Constructors  & Destructors
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
CFD2D::CFD2D () {}
CFD2D::CFD2D (char hoge[64]) {
  sprintf (piyo, "inp.2d.%s", hoge);
  if ((fp_in = fopen (piyo, "r")) == NULL) {
    printf ("No such file as %s(>_<;)\n", piyo);
    exit (1);
  }
  fscanf (fp_in, "%d", &LLX);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &Lx);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%d", &LLY);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &Ly);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &Rho0);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &U0);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &V0);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &P0);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &Gamma); fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &Co);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &Epsilon);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%lf", &RateOblique);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%d", &output_freq);  fgets (mud, 1024, fp_in);
  fscanf (fp_in, "%d", &output_freq_1d);  fgets (mud, 1024, fp_in);

  dx = Lx / (double)LLX;
  dy = Ly / (double)LLY;
  dxi = sqrt (2.0) * dx;
  deta = sqrt (2.0) * dy;
  kappa = 1.0/3.0;
  b = (3.0 - kappa) / (1.0 - kappa);

  n_pulse = 0;


  //  RateOblique = 0.0;
  RateNormal = 1.0 - RateOblique;
  RI2 = 1.0 / sqrt (2.0);

  printf ("##%s loaded\n", piyo);
  printf ("Size of x-stencil:   %d (-)\n", LLX);
  printf ("Length of x-stencil: %lf (m)\n", Lx);
  printf ("Size of x-stencil:   %d (-)\n", LLY);
  printf ("Length of x-stencil: %lf (m)\n", Ly);
  printf ("Rho0:              %lf (kg/m3)\n", Rho0);
  printf ("U0:                %lf (m/s)\n", U0);
  printf ("P0:                %lf (Pa)\n", P0);
  printf ("Gamma:             %lf (-)\n", Gamma);
  printf ("Courant Number:    %lf (-)\n", Co);
  printf ("Epsilon:           %lf (-)\n", Epsilon);
  printf ("dx:                %le (m)\n", dx);
  printf ("kappa              %le (-)\n", kappa);
  printf ("b:                 %le (-)\n", b);
  printf ("normal co-ordinate:  %lf\n", RateNormal);
  printf ("oblique co-ordinate:  %lf\n", RateOblique);
  printf ("output_frequcncy:  %d\n", output_freq);
  printf ("1d output_frequcncy:  %d\n\n", output_freq_1d);
  fclose (fp_in);

  Rho = new double*[LLX];
  U = new double*[LLX];
  V = new double*[LLX];
  P = new double*[LLX];
  RhoU = new double*[LLX];
  RhoV = new double*[LLX];
  RhoE = new double*[LLX];
  H = new double*[LLX];
  C = new double*[LLX];
  RhoL = new double*[LLX];
  RhoR = new double*[LLX];
  UL = new double*[LLX];
  UR = new double*[LLX];
  VL = new double*[LLX];
  VR = new double*[LLX];
  HL = new double*[LLX];
  HR = new double*[LLX];
  PL = new double*[LLX];
  PR = new double*[LLX];
  RhoUL = new double*[LLX];
  RhoUR = new double*[LLX];
  RhoVL = new double*[LLX];
  RhoVR = new double*[LLX];
  RhoEL = new double*[LLX];
  RhoER = new double*[LLX];
  DeltaPlus = new double*[LLX];
  DeltaMinus = new double*[LLX];
  RhoFluxE = new double*[LLX];
  RhoUFluxE = new double*[LLX];
  RhoVFluxE = new double*[LLX];
  RhoEFluxE = new double*[LLX];
  RhoFluxF = new double*[LLX];
  RhoUFluxF = new double*[LLX];
  RhoVFluxF = new double*[LLX];
  RhoEFluxF = new double*[LLX];
  RoeRho = new double*[LLX];
  RoeU = new double*[LLX];
  RoeV = new double*[LLX];
  RoeH = new double*[LLX];
  RoeC = new double*[LLX];
  Rhotmp = new double*[LLX];
  RhoUtmp = new double*[LLX];
  RhoVtmp = new double*[LLX];
  RhoEtmp = new double*[LLX];
  MuE = new double*[LLX];
  RhoS = new double*[LLX];
  RhoUS = new double*[LLX];
  RhoVS = new double*[LLX];
  RhoES = new double*[LLX];
  TauXX = new double*[LLX];
  TauXY = new double*[LLX];
  TauYY = new double*[LLX];
  SigmaXX = new double*[LLX];
  SigmaYY = new double*[LLX];
  Mask = new double *[LLX];
  MaskSurface = new double *[LLX];
  MaskNormalX = new double *[LLX];
  MaskNormalY = new double *[LLX];
  RhoFB = new double *[LLX];
  RhoUFB = new double *[LLX];
  RhoVFB = new double *[LLX];
  RhoEFB = new double *[LLX];
  for (i = 0; i < LLX; i++) {
    Rho[i] = new double[LLY];
    U[i] = new double[LLY];
    V[i] = new double[LLY];
    P[i] = new double[LLY];
    RhoU[i] = new double[LLY];
    RhoV[i] = new double[LLY];
    RhoE[i] = new double[LLY];
    H[i] = new double[LLY];
    C[i] = new double[LLY];
    RhoL[i] = new double[LLY];
    RhoR[i] = new double[LLY];
    UL[i] = new double[LLY];
    UR[i] = new double[LLY];
    VL[i] = new double[LLY];
    VR[i] = new double[LLY];
    HL[i] = new double[LLY];
    HR[i] = new double[LLY];
    PR[i] = new double[LLY];
    PL[i] = new double[LLY];
    RhoUL[i] = new double[LLY];
    RhoUR[i] = new double[LLY];
    RhoVL[i] = new double[LLY];
    RhoVR[i] = new double[LLY];
    RhoEL[i] = new double[LLY];
    RhoER[i] = new double[LLY];
    DeltaPlus[i] = new double[LLY];
    DeltaMinus[i] = new double[LLY];
    RhoFluxE[i] = new double[LLY];
    RhoUFluxE[i] = new double[LLY];
    RhoVFluxE[i] = new double[LLY];
    RhoEFluxE[i] = new double[LLY];
    RhoFluxF[i] = new double[LLY];
    RhoUFluxF[i] = new double[LLY];
    RhoVFluxF[i] = new double[LLY];
    RhoEFluxF[i] = new double[LLY];
    RoeRho[i] = new double[LLY];
    RoeU[i] = new double[LLY];
    RoeV[i] = new double[LLY];
    RoeH[i] = new double[LLY];
    RoeC[i] = new double[LLY];
    Rhotmp[i] = new double[LLY];
    RhoUtmp[i] = new double[LLY];
    RhoVtmp[i] = new double[LLY];
    RhoEtmp[i] = new double[LLY];
    MuE[i] = new double[LLY];
    RhoS[i] = new double[LLY];
    RhoUS[i] = new double[LLY];
    RhoVS[i] = new double[LLY];
    RhoES[i] = new double[LLY];
    TauXX[i] = new double[LLY];
    TauXY[i] = new double[LLY];
    TauYY[i] = new double[LLY];
    SigmaXX[i] = new double[LLY];
    SigmaYY[i] = new double[LLY];
    Mask[i] = new double[LLY];
    MaskSurface[i] = new double[LLY];
    MaskNormalX[i] = new double [LLY];
    MaskNormalY[i] = new double [LLY];
    RhoFB[i] = new double[LLY];
    RhoUFB[i] = new double[LLY];
    RhoVFB[i] = new double[LLY];
    RhoEFB[i] = new double[LLY];
  }
  //  sprintf (piyo, "%s", hoge);


  //  test for mask
  //  fp_out = fopen ("test_mask.dat", "w");
  fp_out = fopen ("inp.mask", "w");
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      if (pow (i - 50, 2) + pow (j - 50, 2) < pow (50.0, 2)) {
	Mask[i][j] = 1.0;
      } else {
	Mask[i][j] = 0.0;
      }
      fprintf (fp_out, "%d\t%d\t%lf\n", i, j, Mask[i][j]);
    }
  }
  fclose (fp_out);

  //  LoadMask ();
}


CFD2D::~CFD2D () {

  for (i = 0; i < LLX; i++) {
    delete [] Rho[i];
    delete [] U[i];
    delete [] V[i];
    delete [] P[i];
    delete [] RhoU[i];
    delete [] RhoV[i];
    delete [] RhoE[i];
    delete [] H[i];
    delete [] C[i];
    delete [] RhoL[i];
    delete [] RhoR[i];
    delete [] UL[i];
    delete [] UR[i];
    delete [] VL[i];
    delete [] VR[i];
    delete [] HL[i];
    delete [] HR[i];
    delete [] PR[i];
    delete [] PL[i];
    delete [] RhoUL[i];
    delete [] RhoUR[i];
    delete [] RhoVL[i];
    delete [] RhoVR[i];
    delete [] RhoEL[i];
    delete [] RhoER[i];
    delete [] DeltaPlus[i];
    delete [] DeltaMinus[i];
    delete [] RhoFluxE[i];
    delete [] RhoUFluxE[i];
    delete [] RhoVFluxE[i];
    delete [] RhoEFluxE[i];
    delete [] RhoFluxF[i];
    delete [] RhoUFluxF[i];
    delete [] RhoVFluxF[i];
    delete [] RhoEFluxF[i];
    delete [] RoeRho[i];
    delete [] RoeU[i];
    delete [] RoeV[i];
    delete [] RoeH[i];
    delete [] RoeC[i];
    delete [] Rhotmp[i];
    delete [] RhoUtmp[i];
    delete [] RhoVtmp[i];
    delete [] RhoEtmp[i];
    delete [] MuE[i];
    delete [] RhoS[i];
    delete [] RhoUS[i];
    delete [] RhoVS[i];
    delete [] RhoES[i];
    delete [] TauXX [i];
    delete [] TauXY [i];
    delete [] TauYY [i];
    delete [] SigmaXX [i];
    delete [] SigmaYY [i];
    delete [] Mask [i];
    delete [] MaskSurface [i];
    delete [] MaskNormalX [i];
    delete [] MaskNormalY [i];
    delete [] RhoFB [i];
    delete [] RhoUFB [i];
    delete [] RhoVFB [i];
    delete [] RhoEFB [i];
  }

  delete [] Rho;
  delete [] U;
  delete [] V;
  delete [] P;
  delete [] RhoU;
  delete [] RhoV;
  delete [] RhoE;
  delete [] H;
  delete [] C;
  delete [] RhoL;
  delete [] RhoR;
  delete [] UL;
  delete [] UR;
  delete [] VL;
  delete [] VR;
  delete [] HL;
  delete [] HR;
  delete [] PL;
  delete [] PR;
  delete [] RhoUL;
  delete [] RhoUR;
  delete [] RhoVL;
  delete [] RhoVR;
  delete [] RhoEL;
  delete [] RhoER;
  delete [] DeltaPlus;
  delete [] DeltaMinus;
  delete [] RhoFluxE;
  delete [] RhoUFluxE;
  delete [] RhoVFluxE;
  delete [] RhoEFluxE;
  delete [] RhoFluxF;
  delete [] RhoUFluxF;
  delete [] RhoVFluxF;
  delete [] RhoEFluxF;
  delete [] RoeRho;
  delete [] RoeU;
  delete [] RoeV;
  delete [] RoeH;
  delete [] RoeC;
  delete [] Rhotmp;
  delete [] RhoUtmp;
  delete [] RhoVtmp;
  delete [] RhoEtmp;
  delete [] MuE;
  delete [] RhoS;
  delete [] RhoUS;
  delete [] RhoVS;
  delete [] RhoES;
  delete [] TauXX;
  delete [] TauXY;
  delete [] TauYY;
  delete [] SigmaXX;
  delete [] SigmaYY;
  delete [] Mask;
  delete [] MaskSurface;
  delete [] MaskNormalX;
  delete [] MaskNormalY;
  delete [] RhoFB;
  delete [] RhoUFB;
  delete [] RhoVFB;
  delete [] RhoEFB;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////    End Constructors & Destructors
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////
//        MUSCL Scheme with Roe Average         //
//////////////////////////////////////////////////
double CFD2D::RoeMUSCL (int n, double dt) {
  //  Unknown Vars.
  for (i = 0; i < LLX; i++)  {
    for (j = 0; j < LLY; j++)
      RhoU[i][j] = Rho[i][j] * U[i][j];
  }
  for (i = 0; i < LLX; i++)  {
    for (j = 0; j < LLY; j++)
      RhoV[i][j] = Rho[i][j] * V[i][j];
  }
  for (i = 0; i < LLX; i++)  {
    for (j = 0; j < LLY; j++) 
      RhoE[i][j] = P[i][j] / (Gamma - 1.0) + 0.5 * Rho[i][j] * (pow (U[i][j], 2) + pow (V[i][j], 2));
  }

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++)
      H[i][j] = Gamma / (Gamma - 1.0) * P[i][j] / Rho[i][j] + 0.5 * (pow (U[i][j], 2) + pow (V[i][j], 2));
  }
  for (i = 0; i < LLX; i++)  {
    for (j = 0; j < LLY; j++) 
      C[i][j] = sqrt (fabs (Gamma * P[i][j] / Rho[i][j]));
    if (isnan (C[i][j])) {
      printf ("Sound Speed nan\n (i, j, Rho, P, Gamma, C) = (%d, %d, %le, %le, %le, %le)\n", i, j, Rho[i][j], P[i][j], Gamma, C[i][j]);
      exit (1);
    }
  }

  //  CFL condition
  CFL ();
  dt = Co * dx / MaxPhaseSpeed;
  cx = dt / dx;
  cy = dt / dy;
  cxi = cx / sqrt (2.0);
  ceta = cy / sqrt (2.0);

  //  printf ("MaxC = %le, dt = %le, cx = %le, cxi = %le\n", MaxC, dt, cx, cxi);
//  printf ("n = %d, PhaseSpeed = %le\n", n, MaxPhaseSpeed);
  if (n % output_freq == 0) {
    printf ("n = %d, PhaseSpeed = %le\n", n, MaxPhaseSpeed);
    //  printf ("dx = %lf, dy = %lf, dt = %lf, Co = %lf\n", dx, dy, dt, Co);
    printf ("dt = %lf, dx = %lf, cx = %lf, cxi = %lf\n", dt, dx, cx, cxi);
  }


  //  Normal direction numerical flux  
  QLeftRightX ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoEL
  RoeAverage ();
  XDirectionalFlux ();
  QLeftRightY ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoEL
  RoeAverage ();
  YDirectionalFlux ();

  printf ("numerical flux for X-Y is obtained\n");

  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      Rhotmp[i][j] = Rho[i][j] - RateNormal * (cx * RhoFluxE[i][j] + cy * RhoFluxF[i][j]); 
      //      Rhotmp[i][j] = Rho[i][j] - (cx * RhoFluxE[i][j] + cy * RhoFluxF[i][j]); 
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      RhoUtmp[i][j] = RhoU[i][j] - RateNormal * (cx * RhoUFluxE[i][j] + cy * RhoUFluxF[i][j]); 
      //      RhoUtmp[i][j] = RhoU[i][j] - (cx * RhoUFluxE[i][j] + cy * RhoUFluxF[i][j]); 
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      RhoVtmp[i][j] = RhoV[i][j] - RateNormal * (cx * RhoVFluxE[i][j] + cy * RhoVFluxF[i][j]); 
      //      RhoVtmp[i][j] = RhoV[i][j] - (cx * RhoVFluxE[i][j] + cy * RhoVFluxF[i][j]); 
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      RhoEtmp[i][j] = RhoE[i][j] - RateNormal * (cx * RhoEFluxE[i][j] + cy * RhoEFluxF[i][j]); 
      //      RhoEtmp[i][j] = RhoE[i][j] - (cx * RhoEFluxE[i][j] + cy * RhoEFluxF[i][j]); 
    }
  }




  //  MultipleDirection Numerical Flux
  QLeftRightXi ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoELn
  RoeAverage ();
  XiDirectionalFlux ();
  QLeftRightEta ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoEL
  RoeAverage ();
  EtaDirectionalFlux ();

  printf ("numerical flux for Xi-Eta is obtained\n");


  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      //      Rhotmp[i][j] = Rho[i][j] - cx * RhoFluxE[i][j] - cy * RhoFluxF[i][j] ; 
      //      Rhotmp[i][j] = Rhotmp[i][j] - RateOblique * (cxi * RhoFluxE[i][j] + ceta * RhoFluxF[i][j]); w
      Rhotmp[i][j] = Rhotmp[i][j] - RI2 * RateOblique * (cxi * RhoFluxE[i][j] + ceta * RhoFluxF[i][j]); 
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      //      RhoUtmp[i][j] = RhoU[i][j] - cx * RhoUFluxE[i][j] - cy * RhoUFluxF[i][j]; 
      RhoUtmp[i][j] = RhoUtmp[i][j] - RI2 * RateOblique * (cxi * RhoUFluxE[i][j] + ceta * RhoUFluxF[i][j]); 
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      //      RhoVtmp[i][j] = RhoV[i][j] - cx * RhoVFluxE[i][j] - cy * RhoVFluxF[i][j]; 
      RhoVtmp[i][j] = RhoVtmp[i][j] - RI2 * RateOblique * (cxi * RhoVFluxE[i][j] + ceta * RhoVFluxF[i][j]); 
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      //      RhoEtmp[i][j] = RhoE[i][j] - cx * RhoEFluxE[i][j] - cy * RhoEFluxF[i][j]; 
      RhoEtmp[i][j] = RhoEtmp[i][j] - RI2 * RateOblique * (cxi * RhoEFluxE[i][j] + ceta * RhoEFluxF[i][j]);  
   }
  }

  printf ("new vars. obtained\n");


  //////////////////////////////////////////////////
  //             call  LES routine                //
  //////////////////////////////////////////////////
  LES (n, dt);


  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      Rho[i][j] = Rhotmp[i][j];
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      RhoU[i][j] = RhoUtmp[i][j] + RhoUS[i][j];
    }
  } 
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      RhoV[i][j] = RhoVtmp[i][j] + RhoVS[i][j];
    }
  }
  for (i = 3; i < LLX-3; i++) {
    for (j = 3; j < LLY-3; j++) {
      RhoE[i][j] = RhoEtmp[i][j] + RhoES[i][j];
    }
  }



  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      U[i][j] = RhoU[i][j] / Rho[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      V[i][j] = RhoV[i][j] / Rho[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      P[i][j] = (RhoE[i][j] - 0.5*Rho[i][j]*(pow(U[i][j],2) + pow(V[i][j], 2))) * (Gamma - 1.0);
    }
  }


  printf ("%d th iteration completed\n", n);
  
  return Co * dx / MaxPhaseSpeed;  
}


//////////////////////////////////////////////////
//                     CFL                      //
//////////////////////////////////////////////////
inline void CFD2D::CFL () {
  MaxPhaseSpeed = 0.0;
  MaxU = 0.0; MaxV = 0.0; MaxC = 0.0;  
  MinP = 0.0;

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) 
      MaxU = fmax (MaxU, fabs(U[i][j]));
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) 
      MaxV = fmax (MaxV, fabs (V[i][j]));
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) 
      MaxC = fmax (MaxC, fabs(C[i][j]));
  }
  MaxPhaseSpeed = MaxU + MaxV + MaxC + 1.0;


  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      MinP = fmin (MinP, fabs(P[i][j]));
      if (MinP < 0.0) {
 	printf ("P to be a negative value at (i, j) = (%d, %d).\n. The application going to be shut down(>_<;)\n", i, j);
 	exit (1);
      }
    }
  }


}


//////////////////////////////////////////////////
//                 X-direction                  //
//////////////////////////////////////////////////
inline void CFD2D::XDirectionalFlux () {
  for (i = 2; i < LLX-2; i++) {
    for (j = 2; j < LLY-2; j++) {
      //      printf ("hoge i, j = %d, %d\n", i, j);
      //      Obtain Flux of X-directional Emptyset|i+1/2
      RightDiagonalMatrixX (i, j);
      LeftDiagonalMatrixX (i, j);
      //      printf ("Diagonal Matrix for X direction obtained\n");

      AbsLambdaX (i, j);
      AbsJacobiAXY ();
      //      printf ("Jacobi matrix obtained X-direction obtained\n");

      RhoOutE = RhoFluxX (i, j);
      RhoUOutE = RhoUFluxX (i, j);
      RhoVOutE = RhoVFluxX (i, j);
      RhoEOutE = RhoEFluxX (i, j);

      //      Obtain Flux of X-directional Emptyset|i-1/2
      //      RoeAverage (i-1, j);
      //      RoeAverage ();
      RightDiagonalMatrixX (i-1, j);
      LeftDiagonalMatrixX (i-1, j);
      AbsLambdaX (i-1, j);
      AbsJacobiAXY ();

      RhoInE = RhoFluxX (i-1, j);
      RhoUInE = RhoUFluxX (i-1, j);
      RhoVInE = RhoVFluxX (i-1, j);
      RhoEInE = RhoEFluxX (i-1, j);


      RhoFluxE[i][j]  = RhoOutE - RhoInE;
      RhoUFluxE[i][j] = RhoUOutE - RhoUInE;
      RhoVFluxE[i][j] = RhoVOutE - RhoVInE;
      RhoEFluxE[i][j] = RhoEOutE - RhoEInE;
    }
  }


}
inline void CFD2D::QLeftRightX () {
   EmptySetX (Rho, DeltaPlus, DeltaMinus, RhoL, RhoR);
   EmptySetX (U, DeltaPlus, DeltaMinus, UL, UR);
   EmptySetX (V, DeltaPlus, DeltaMinus, VL, VR);
   EmptySetX (P, DeltaPlus, DeltaMinus, PL, PR);
   EmptySetX (H, DeltaPlus, DeltaMinus, HL, HR);

   //  Rho|i+1/2, U|i+1/2, V|i+1/2 --> RhoUR, RhoUL, RhoVR, RhoVL
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoUL[i][j] = RhoL[i][j] * UL[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoUR[i][j] = RhoR[i][j] * UR[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoVL[i][j] = RhoL[i][j] * VL[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoVR[i][j] = RhoR[i][j] * VR[i][j];
     }
   }

  //   Rho|i+1/2, P|i+1/2, U|i+1/2, V|i+1/2 -->
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoEL[i][j] = PL[i][j] / (Gamma - 1.0) + 0.5 * RhoL[i][j] * (pow (UL[i][j], 2) + pow (VL[i][j], 2));
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoER[i][j] = PR[i][j] / (Gamma - 1.0) + 0.5 * RhoR[i][j] * (pow (UR[i][j], 2) + pow (VR[i][j], 2));
    }
  }
}



//void CFD2D::EmptySetLRX (double** Plus, double** DeltaMinus, double** EmptySet, double** EmptySetL, double** EmptySetR) {
//void CFD2D::EmptySetX (double** In, double** Plus, double** Minus, double** EmptySet, double** EmptySetL, double** EmptySetR) {
inline void CFD2D::EmptySetX (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR) {
  for (i = 1; i < LLX-1; i++) {
    for (j = 0; j < LLY; j++) {
      Plus[i][j]  = minmod (EmptySet[i+1][j] - EmptySet[i][j], b*(EmptySet[i][j] - EmptySet[i-1][j]));

    }
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 0; j < LLY; j++) {
      Minus[i][j] = minmod (EmptySet[i][j] - EmptySet[i-1][j], b*(EmptySet[i+1][j] - EmptySet[i][j]));
      //      printf ("%d, %d, Rho = %le, DeltaMinus = %le\n", 
      //	      i, j, In[i][j], DeltaMinus[i][j]);
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 0; j < LLY; j++) {
      s = (2.0 * Plus[i][j] * Minus[i][j] + Epsilon) / (pow (Plus[i][j], 2) + pow (Minus[i][j], 2) + Epsilon);
      EmptySetL[i][j] = EmptySet[i][j] + 0.25*s * ((1-kappa*s)* Plus[i][j] + (1+kappa*s) * Minus[i][j]); 
      //      printf ("%d, %d, Rho = %le, RhoL = %le\n", 
      //      	      i, j, EmptySet[i][j], EmptySetL[i][j]);
    }
  }  
  for (i = 1; i < LLX-1; i++) {
    for (j = 0; j < LLY; j++) {
      s = (2.0 * Plus[i+1][j] * Minus[i+1][j] + Epsilon) / (pow (Plus[i+1][j], 2) + pow (Minus[i+1][j], 2) + Epsilon);
      EmptySetR[i][j] = EmptySet[i+1][j] - 0.25*s * ((1-kappa*s) * Plus[i+1][j] + (1+kappa*s) * Minus[i+1][j]); 
    }
  }


}

inline double CFD2D::RhoFluxX (int i, int j) {
  return 0.5 * (RhoUR[i][j] + RhoUL[i][j]
 		       - (AbsJacobiA[0][0] * (RhoR[i][j] - RhoL[i][j]) 
 			  + AbsJacobiA[0][1] * (RhoUR[i][j] - RhoUL[i][j]) 
 			  + AbsJacobiA[0][2] * (RhoVR[i][j] - RhoVL[i][j]) 
 			  + AbsJacobiA[0][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoUFluxX (int i, int j) {
  return 0.5 * (RhoR[i][j] * pow (UR[i][j], 2) + PR[i][j] 
			+ RhoL[i][j] * pow (UL[i][j], 2) + PL[i][j]
			- (AbsJacobiA[1][0] * (RhoR[i][j] - RhoL[i][j]) 
			   + AbsJacobiA[1][1] * (RhoUR[i][j] - RhoUL[i][j]) 
			   + AbsJacobiA[1][2] * (RhoVR[i][j] - RhoVL[i][j]) 
			   + AbsJacobiA[1][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoVFluxX (int i, int j) {
  return 0.5 * (RhoR[i][j] * UR[i][j] * VR[i][j] + RhoL[i][j] * UL[i][j] * VL[i][j]
		       - (AbsJacobiA[2][0] * (RhoR[i][j] - RhoL[i][j])
			  + AbsJacobiA[2][1] * (RhoUR[i][j] - RhoUL[i][j])
			  + AbsJacobiA[2][2] * (RhoVR[i][j] - RhoVL[i][j])
			  + AbsJacobiA[2][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}
inline double CFD2D::RhoEFluxX (int i, int j) {
  return 0.5 * (RhoUR[i][j] * HR[i][j]
		+ RhoUL[i][j] * HL[i][j]
		- (AbsJacobiA[3][0] * (RhoR[i][j] - RhoL[i][j])
		   + AbsJacobiA[3][1] * (RhoUR[i][j] - RhoUL[i][j])
		   + AbsJacobiA[3][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[3][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}




//////////////////////////////////////////////////
//               Roe average                    //
//////////////////////////////////////////////////
inline void CFD2D::RoeAverage () {
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      RoeRho[i][j] = sqrt (RhoR[i][j] * RhoL[i][j]);
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      RoeU[i][j]   = (sqrt(RhoR[i][j])*UR[i][j] + sqrt(RhoL[i][j])*UL[i][j]) / (sqrt(RhoL[i][j]) + sqrt(RhoR[i][j]));
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      RoeV[i][j]   = (sqrt(RhoR[i][j])*VR[i][j] + sqrt(RhoL[i][j])*VL[i][j]) / (sqrt(RhoL[i][j]) + sqrt(RhoR[i][j]));
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      RoeH[i][j]   = (sqrt(RhoR[i][j])*HR[i][j] + sqrt(RhoL[i][j])*HL[i][j]) / (sqrt(RhoL[i][j]) + sqrt(RhoR[i][j]));
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      RoeC[i][j]   = sqrt(fabs((Gamma - 1.0) * (RoeH[i][j] - 0.5 * pow (RoeU[i][j], 2))));
    }
  }
    //    printf ("RoeRho = %lf, RoeU = %lf, RoeH = %lf, RoeC = %lf\n", RoeRho, RoeU, RoeH, RoeC);
}

//////////////////////////////////////////////////
//     X-directional Jacobi Matrix et al.       //
//////////////////////////////////////////////////
inline void CFD2D::RightDiagonalMatrixX (int i, int j) {
//     int ii, kk; 

  //    for (j = 0; j < 3; j++)
  //      DiagR[j][0] = 1.0;
  DiagR[0][0] = 1.0;
  DiagR[0][1] = 1.0;
  DiagR[0][2] = 1.0;
  DiagR[0][3] = 0.0;

  DiagR[1][0] = RoeU[i][j] - RoeC[i][j];   
  DiagR[1][1] = RoeU[i][j];   
  DiagR[1][2] = RoeU[i][j] + RoeC[i][j];  
  DiagR[1][3] = 0.0;

  DiagR[2][0] = RoeV[i][j];    
  DiagR[2][1] = RoeV[i][j];    
  DiagR[2][2] = RoeV[i][j];    
  DiagR[2][3] = 1.0;  

  DiagR[3][0] = RoeH[i][j] - RoeC[i][j] * RoeU[i][j];
  DiagR[3][1] = 0.5 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));
  DiagR[3][2] = RoeH[i][j] + RoeC[i][j] * RoeU[i][j];
  DiagR[3][3] = RoeV[i][j];    

//    for (ii = 0; ii < 3; ii++) {
//      for (kk = 0; kk < 3; kk++)  printf ("%lf\t ", DiagR[ii][kk]);
//      printf ("\n");
//    }
//     printf ("\n");

}



inline void CFD2D::LeftDiagonalMatrixX (int i, int j) {
    //  Compose Left Diagonal Matrix for F|i-1/2
  //  Eta2 = Eta1 * (pow (RoeU, 2) + pow (RoeV, 2));
  Eta2 = (Gamma - 1.0) / pow (RoeC[i][j], 2);
  Eta1 = 0.5 * Eta2 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));

  DiagL[0][0] = 0.5 * (Eta1 + RoeU[i][j] / RoeC[i][j]);
  DiagL[0][1] = -0.5 * (Eta2 * RoeU[i][j] + 1.0 / RoeC[i][j]);
  DiagL[0][2] = -0.5 * Eta2 * RoeV[i][j];
  DiagL[0][3] = 0.5 * Eta2;

  DiagL[1][0] = 1.0 - Eta1;
  DiagL[1][1] = Eta2 * RoeU[i][j];
  DiagL[1][2] = Eta2 * RoeV[i][j];
  DiagL[1][3] = -Eta2;

  DiagL[2][0] = 0.5 * (Eta1 - RoeU[i][j] / RoeC[i][j]);
  DiagL[2][1] = 0.5 * (1.0 / RoeC[i][j] - Eta2 * RoeU[i][j]);
  DiagL[2][2] = -0.5 * Eta2 * RoeV[i][j];
  DiagL[2][3] = 0.5 * Eta2;

  DiagL[3][0] = -RoeV[i][j];
  DiagL[3][1] = 0.0;
  DiagL[3][2] = 1.0;
  DiagL[3][3] = 0.0;
}




inline void CFD2D::AbsLambdaX (int i, int j) {
//     Lambda[0][0] = (RoeU - RoeC);   Lambda[0][1] = 0.0;           Lambda[0][2] = 0.0;                Lambda[0][3] = 0.0;
//     Lambda[1][0] = 0.0;             Lambda[1][1] = fabs (RoeU);   Lambda[1][2] = 0.0;                Lambda[1][3] = 0.0;
//     Lambda[2][0] = 0.0;                  Lambda[2][1] = 0.0;           Lambda[2][2] = fabs (RoeU + RoeC); Lambda[2][3] = 0.0;
//    Lambda[3][0] = 0.0;                  Lambda[3][1] = 0.0;           Lambda[3][2] = 0.0;                Lambda[3][3] = fabs (RoeU);

  AbsLambda[0][0] = fabs(RoeU[i][j] - RoeC[i][j]); AbsLambda[0][1] = 0.0;               AbsLambda[0][2] = 0.0;                             AbsLambda[0][3] = 0.0;
  AbsLambda[1][0] = 0.0;                           AbsLambda[1][1] = fabs (RoeU[i][j]); AbsLambda[1][2] = 0.0;                             AbsLambda[1][3] = 0.0;
  AbsLambda[2][0] = 0.0;                           AbsLambda[2][1] = 0.0;               AbsLambda[2][2] = fabs (RoeU[i][j] + RoeC[i][j]);  AbsLambda[2][3] = 0.0;
  AbsLambda[3][0] = 0.0;                           AbsLambda[3][1] = 0.0;               AbsLambda[3][2] = 0.0;                             AbsLambda[3][3] = fabs (RoeU[i][j]);

//     AbsLambda[0][0] = (RoeU[i][j] - RoeC[i][j]); AbsLambda[0][1] = 0.0;          AbsLambda[0][2] = 0.0;                        AbsLambda[0][3] = 0.0;
//     AbsLambda[1][0] = 0.0;                       AbsLambda[1][1] = (RoeU[i][j]); AbsLambda[1][2] = 0.0;                        AbsLambda[1][3] = 0.0;
//     AbsLambda[2][0] = 0.0;                       AbsLambda[2][1] = 0.0;          AbsLambda[2][2] = (RoeU[i][j] + RoeC[i][j]);  AbsLambda[2][3] = 0.0;
//     AbsLambda[3][0] = 0.0;                       AbsLambda[3][1] = 0.0;          AbsLambda[3][2] = 0.0;                        AbsLambda[3][3] = (RoeU[i][j]);

}





inline double CFD2D::RhoFluxY (int i, int j) {
  return 0.5 * (RhoVR[i][j] + RhoVL[i][j]
 		       - (AbsJacobiA[0][0] * (RhoR[i][j] - RhoL[i][j]) 
 			  + AbsJacobiA[0][1] * (RhoUR[i][j] - RhoUL[i][j]) 
 			  + AbsJacobiA[0][2] * (RhoVR[i][j] - RhoVL[i][j]) 
 			  + AbsJacobiA[0][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoUFluxY (int i, int j) {
  return 0.5 * (RhoR[i][j] * UR[i][j] * VR[i][j] + RhoL[i][j] * UL[i][j] * VL[i][j]
 			- (AbsJacobiA[1][0] * (RhoR[i][j] - RhoL[i][j]) 
 			   + AbsJacobiA[1][1] * (RhoUR[i][j] - RhoUL[i][j]) 
 			   + AbsJacobiA[1][2] * (RhoVR[i][j] - RhoVL[i][j]) 
 			   + AbsJacobiA[1][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoVFluxY (int i, int j) {
return 0.5 * (RhoR[i][j] * pow (VR[i][j], 2) + PR[i][j] 
 			 + RhoL[i][j] * pow (VL[i][j], 2) + PL[i][j]
 			-(AbsJacobiA[2][0] * (RhoR[i][j] - RhoL[i][j]) 
 			  + AbsJacobiA[2][1] * (RhoUR[i][j] - RhoUL[i][j]) 
 			  + AbsJacobiA[2][2] * (RhoVR[i][j] - RhoVL[i][j]) 
 			  + AbsJacobiA[2][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoEFluxY (int i, int j) {
  return 0.5 * (RhoVR[i][j] * HR[i][j]
		+ RhoVL[i][j] * HL[i][j]
		- (AbsJacobiA[3][0] * (RhoR[i][j] - RhoL[i][j])
		   + AbsJacobiA[3][1] * (RhoUR[i][j] - RhoUL[i][j])
		   + AbsJacobiA[3][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[3][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}

inline void CFD2D::AbsJacobiAXY () {
  //    int j, k, l;
  //  int j, k, l;

  //  Compose Jacobi Matrix
  for (jm = 0; jm < 4; jm++) {
    for (km = 0; km < 4; km++) AbsJacobiAtmp[jm][km] = 0.0;
  }
  for (jm = 0; jm < 4; jm++) {
    for (km = 0; km < 4; km++) AbsJacobiA[jm][km] = 0.0;
  }
  for (jm = 0; jm < 4; jm++) {
    for (km = 0; km < 4; km++) {
      for (lm = 0; lm < 4; lm++) AbsJacobiAtmp[jm][km] += DiagR[jm][lm] * AbsLambda[lm][km];
    }
  }
  for (jm = 0; jm < 4; jm++) {
    for (km = 0; km < 4; km++) {
      for (lm = 0; lm < 4; lm++) AbsJacobiA[jm][km] += AbsJacobiAtmp[jm][lm] * DiagL[lm][km];
      //      if (isnan (AbsJacobiA[jm][km])) printf ("(i, j, jm, km) = (%d, %d, %d, %d)\n(L, R, Lambda, J) = (%le, %le, %le, %le)\n", 
      //      					      i, j, jm, km, DiagL[jm][km], DiagR[jm][km], AbsLambda[jm][km], AbsJacobiAtmp[jm][km]);
    }
  }

}

//////////////////////////////////////////////////
//                Y-direction                   //
//////////////////////////////////////////////////
inline void CFD2D::YDirectionalFlux () {
  //  QLeftRightY ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoEL
  //  printf ("Emptyset for Y-dire// // // // ction obtained\n");

  //  RoeAverage ();
  //  printf ("Roe average for Y-direction obtained\n");

  for (i = 2; i < LLX-2; i++) {
    for (j = 2; j < LLY-2; j++) {
      //      Obtain of Y-directional Flux Emptyset|j+1/2
      RightDiagonalMatrixY (i, j);
      LeftDiagonalMatrixY (i, j);
      AbsLambdaY (i, j);
      AbsJacobiAXY ();

      RhoOutF = RhoFluxY (i, j);
      RhoUOutF = RhoUFluxY (i, j);
      RhoVOutF = RhoVFluxY (i, j);
      RhoEOutF = RhoEFluxY (i, j);

      //      Obtain of Y-directional Flux Emptyset|j-1/2
      RightDiagonalMatrixY (i, j-1);
      LeftDiagonalMatrixY (i, j-1);
      AbsLambdaY (i, j-1);
      AbsJacobiAXY ();

      RhoInF = RhoFluxY (i, j-1);
      RhoUInF = RhoUFluxY (i, j-1);
      RhoVInF = RhoVFluxY (i, j-1);
      RhoEInF = RhoEFluxY (i, j-1);

      RhoFluxF[i][j]  = RhoOutF - RhoInF;
      RhoUFluxF[i][j] = RhoUOutF - RhoUInF;
      RhoVFluxF[i][j] = RhoVOutF - RhoVInF;
      RhoEFluxF[i][j] = RhoEOutF - RhoEInF;


      //      printf ("RhoF = %lf, RhoUF = %lf, RhoEF = %lf\n", RhoF[i][j], RhoUF[i][j], RhoEF[i][j]);
      //      printf ("RhoER = %lf, RhoEL = %lf\n", RhoER[i][j], RhoEL[i][j]);
      //      printf ("A11 = %le, A12 = %le, A13=%le\n", JacobiA[0][0], JacobiA[0][1], JacobiA[0][2]);
      //      printf ("y, %d, %d, %le, %le\n", i, j, RhoOutF, RhoInF);
    }
  }


}

inline void CFD2D::QLeftRightY () {
  EmptySetY (Rho, DeltaPlus, DeltaMinus, RhoL, RhoR);
  EmptySetY (U, DeltaPlus, DeltaMinus, UL, UR);
  EmptySetY (V, DeltaPlus, DeltaMinus, VL, VR);
  EmptySetY (P, DeltaPlus, DeltaMinus, PL, PR);
  EmptySetY (H, DeltaPlus, DeltaMinus, HL, HR);

  //  RhoR, RhoL, UR, UL, VR, VL --> RhoUR, RhoUL, RhoVL, RhoVR
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      RhoUL[i][j] = RhoL[i][j] * UL[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      RhoUR[i][j] = RhoR[i][j] * UR[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      RhoVL[i][j] = RhoL[i][j] * VL[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      RhoVR[i][j] = RhoR[i][j] * VR[i][j];
    }
  }


  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      RhoEL[i][j] = PL[i][j] / (Gamma - 1.0) + 0.5 * RhoL[i][j] * (pow (UL[i][j], 2) + pow (VL[i][j], 2));
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      RhoER[i][j] = PR[i][j] / (Gamma - 1.0) + 0.5 * RhoR[i][j] * (pow (UR[i][j], 2) + pow (VR[i][j], 2));
    }
  }
}



void CFD2D::EmptySetLRY (double** DeltaPlus, double** DeltaMinus, double** EmptySet, double** EmptySetL, double** EmptySetR) {
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      //      s = (2.0 * DeltaPlus[i][j] * DeltaMinus[i][j] + Epsilon) / (pow (DeltaPlus[i][j], 2) + pow (DeltaMinus[i][j], 2) + Epsilon);
      s = Slope (DeltaPlus[i][j], DeltaMinus[i][j], Epsilon); 
      EmptySetL[i][j] = EmptySet[i][j] + 0.25*s * ((1-kappa*s)* DeltaPlus[i][j] + (1+kappa*s) * DeltaMinus[i][j]); 
    }
  }  
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      //      s = (2.0 * DeltaPlus[i][j+1] * DeltaMinus[i][j+1] + Epsilon) / (pow (DeltaPlus[i][j+1], 2) + pow (DeltaMinus[i][j+1], 2) + Epsilon);
      s = Slope (DeltaPlus[i][j+1], DeltaMinus[i][j+1], Epsilon); 
      EmptySetR[i][j] = EmptySet[i][j+1] - 0.25*s * ((1-kappa*s) * DeltaPlus[i][j+1] + (1+kappa*s) * DeltaMinus[i][j+1]); 
    }
  }

}


void CFD2D::EmptySetY (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR) {

  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) 
      Plus[i][j]  = minmod (EmptySet[i][j+1] - EmptySet[i][j], b*(EmptySet[i][j] - EmptySet[i][j-1]));
  }
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) 
      Minus[i][j] = minmod (EmptySet[i][j] - EmptySet[i][j-1], b*(EmptySet[i][j+1] - EmptySet[i][j]));
  }
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      //      s = (2.0 * DeltaPlus[i][j] * DeltaMinus[i][j] + Epsilon) / (pow (DeltaPlus[i][j], 2) + pow (DeltaMinus[i][j], 2) + Epsilon);
      s = Slope (Plus[i][j], Minus[i][j], Epsilon); 
      EmptySetL[i][j] = EmptySet[i][j] + 0.25*s * ((1-kappa*s)* Plus[i][j] + (1+kappa*s) * Minus[i][j]); 
    }
  }  
  for (i = 0; i < LLX; i++) {
    for (j = 1; j < LLY-1; j++) {
      //      s = (2.0 * DeltaPlus[i][j+1] * DeltaMinus[i][j+1] + Epsilon) / (pow (DeltaPlus[i][j+1], 2) + pow (DeltaMinus[i][j+1], 2) + Epsilon);
      s = Slope (Plus[i][j+1], Minus[i][j+1], Epsilon); 
      EmptySetR[i][j] = EmptySet[i][j+1] - 0.25*s * ((1-kappa*s) * Plus[i][j+1] + (1+kappa*s) * Minus[i][j+1]); 
    }
  }

}



inline void CFD2D::RightDiagonalMatrixY (int i, int j) {
//     int ii, kk; 

  //    for (j = 0; j < 3; j++)
  //      DiagR[j][0] = 1.0;
  DiagR[0][0] = 1.0;
  DiagR[0][1] = 1.0;
  DiagR[0][2] = 1.0;
  DiagR[0][3] = 0.0;

  DiagR[1][0] = RoeU[i][j];    
  DiagR[1][1] = RoeU[i][j];    
  DiagR[1][2] = RoeU[i][j];    
  DiagR[1][3] = 1.0;  

  DiagR[2][0] = RoeV[i][j] - RoeC[i][j];   
  DiagR[2][1] = RoeV[i][j];   
  DiagR[2][2] = RoeV[i][j] + RoeC[i][j];  
  DiagR[2][3] = 0.0;

  DiagR[3][0] = RoeH[i][j] - RoeC[i][j] * RoeV[i][j];
  DiagR[3][1] = 0.5 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));
  DiagR[3][2] = RoeH[i][j] + RoeC[i][j] * RoeV[i][j];
  DiagR[3][3] = RoeU[i][j];    

//    for (ii = 0; ii < 3; ii++) {
//      for (kk = 0; kk < 3; kk++)  printf ("%lf\t ", DiagR[ii][kk]);
//      printf ("\n");
//    }
//     printf ("\n");

}



inline void CFD2D::LeftDiagonalMatrixY (int i, int j) {
    //  Compose Left Diagonal Matrix for F|i-1/2
  //  Eta2 = Eta1 * (pow (RoeU, 2) + pow (RoeV, 2));
  Eta2 = (Gamma - 1.0) / pow (RoeC[i][j], 2);
  Eta1 = 0.5 * Eta2 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));
  

  DiagL[0][0] = 0.5 * (Eta1 + RoeV[i][j] / RoeC[i][j]);
  DiagL[0][1] = -0.5 * Eta2 * RoeU[i][j];
  DiagL[0][2] = -0.5 * (Eta2 * RoeV[i][j] + 1.0 / RoeC[i][j]);
  DiagL[0][3] = 0.5 * Eta2;

  DiagL[1][0] = 1.0 - Eta1;
  DiagL[1][1] = Eta2 * RoeU[i][j];
  DiagL[1][2] = Eta2 * RoeV[i][j];
  DiagL[1][3] = -Eta2;

  DiagL[2][0] = 0.5 * (Eta1 - RoeV[i][j] / RoeC[i][j]);
  DiagL[2][1] = -0.5 * Eta2 * RoeU[i][j];
  DiagL[2][2] = 0.5 * (1.0 / RoeC[i][j] - Eta2 * RoeV[i][j]);
  DiagL[2][3] = 0.5 * Eta2;

  DiagL[3][0] = -RoeU[i][j];
  DiagL[3][1] = 1.0;
  DiagL[3][2] = 0.0;
  DiagL[3][3] = 0.0;

  //  if (isnan (DiagL[0][0])) 
  //    printf ("(i, j, RoeC, Eta1, Eta2) = (%d, %d, %le, %le, %le)\n",
  //	    i, j, RoeC[i][j], Eta1, Eta2);
}






inline void CFD2D::AbsLambdaY (int i, int j) {
  AbsLambda[0][0] = fabs(RoeV[i][j] - RoeC[i][j]); AbsLambda[0][1] = 0.0;               AbsLambda[0][2] = 0.0;                             AbsLambda[0][3] = 0.0;
  AbsLambda[1][0] = 0.0;                           AbsLambda[1][1] = fabs (RoeV[i][j]); AbsLambda[1][2] = 0.0;                             AbsLambda[1][3] = 0.0;
  AbsLambda[2][0] = 0.0;                           AbsLambda[2][1] = 0.0;               AbsLambda[2][2] = fabs (RoeV[i][j] + RoeC[i][j]);  AbsLambda[2][3] = 0.0;
  AbsLambda[3][0] = 0.0;                           AbsLambda[3][1] = 0.0;               AbsLambda[3][2] = 0.0;                             AbsLambda[3][3] = fabs (RoeV[i][j]);
}




inline double CFD2D::Slope (double DeltaPlus, double DeltaMinus, double Epsilon) {
  return (2.0 * DeltaPlus * DeltaMinus + Epsilon) / (pow (DeltaPlus, 2) + pow (DeltaMinus, 2) + Epsilon);
}


//////////////////////////////////////////////////
//        Multiple Directional Differential     //
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//              Multiple Xi-direction            //
//////////////////////////////////////////////////
inline void CFD2D::XiDirectionalFlux () {
  for (i = 2; i < LLX-2; i++) {
    for (j = 2; j < LLY-2; j++) {
      //      QLeftRightXi ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoELn

      //      RoeAverage ();

      //      Obtain Flux of X-directional Emptyset|i+1/2
      RightDiagonalMatrixXi (i, j);
      LeftDiagonalMatrixXi (i, j);
      //      printf ("Diagonal Matrix for X direction obtained\n");

      AbsLambdaX (i, j);
      AbsJacobiAXY ();
      //      printf ("Jacobi matrix obtained X-direction obtained\n");

      RhoOutE = RhoFluxXi (i, j);
      RhoUOutE = RhoUFluxXi (i, j);
      RhoVOutE = RhoVFluxXi (i, j);
      RhoEOutE = RhoEFluxXi (i, j);

      //      Obtain Flux of X-directional Emptyset|i-1/2
      //      RoeAverage (i-1, j);
      //      RoeAverage ();
      RightDiagonalMatrixX (i-1, j-1);
      LeftDiagonalMatrixX (i-1, j-1);
      AbsLambdaX (i-1, j-1);
      AbsJacobiAXY ();

      RhoInE = RhoFluxXi (i-1, j-1);
      RhoUInE = RhoUFluxXi (i-1, j-1);
      RhoVInE = RhoVFluxXi (i-1, j-1);
      RhoEInE = RhoEFluxXi (i-1, j-1);


      RhoFluxE[i][j]  = RhoOutE - RhoInE;
      RhoUFluxE[i][j] = RhoUOutE - RhoUInE;
      RhoVFluxE[i][j] = RhoVOutE - RhoVInE;
      RhoEFluxE[i][j] = RhoEOutE - RhoEInE;
    }
  }
}
inline void CFD2D::QLeftRightXi () {
   EmptySetXi (Rho, DeltaPlus, DeltaMinus, RhoL, RhoR);
   EmptySetXi (U, DeltaPlus, DeltaMinus, UL, UR);
   EmptySetXi (V, DeltaPlus, DeltaMinus, VL, VR);
   EmptySetXi (P, DeltaPlus, DeltaMinus, PL, PR);
   EmptySetXi (H, DeltaPlus, DeltaMinus, HL, HR);

   //  Rho|i+1/2, U|i+1/2, V|i+1/2 --> RhoUR, RhoUL, RhoVR, RhoVL
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoUL[i][j] = RhoL[i][j] * UL[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoUR[i][j] = RhoR[i][j] * UR[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoVL[i][j] = RhoL[i][j] * VL[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoVR[i][j] = RhoR[i][j] * VR[i][j];
     }
   }

  //   Rho|i+1/2, P|i+1/2, U|i+1/2, V|i+1/2 -->
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoEL[i][j] = PL[i][j] / (Gamma - 1.0) + 0.5 * RhoL[i][j] * (pow (UL[i][j], 2) + pow (VL[i][j], 2));
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoER[i][j] = PR[i][j] / (Gamma - 1.0) + 0.5 * RhoR[i][j] * (pow (UR[i][j], 2) + pow (VR[i][j], 2));
    }
  }
}

inline void CFD2D::EmptySetXi (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR) {
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      Plus[i][j]  = minmod (EmptySet[i+1][j+1] - EmptySet[i][j], b*(EmptySet[i][j] - EmptySet[i-1][j-1]));
    }
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      Minus[i][j] = minmod (EmptySet[i][j] - EmptySet[i-1][j-1], b*(EmptySet[i+1][j+1] - EmptySet[i][j]));
      //      printf ("%d, %d, Rho = %le, DeltaMinus = %le\n", 
      //	      i, j, In[i][j], DeltaMinus[i][j]);
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      s = (2.0 * Plus[i][j] * Minus[i][j] + Epsilon) / (pow (Plus[i][j], 2) + pow (Minus[i][j], 2) + Epsilon);
      EmptySetL[i][j] = EmptySet[i][j] + 0.25*s * ((1-kappa*s)* Plus[i][j] + (1+kappa*s) * Minus[i][j]); 
      //      printf ("%d, %d, Rho = %le, RhoL = %le\n", 
      //      	      i, j, EmptySet[i][j], EmptySetL[i][j]);
    }
  }  
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      s = (2.0 * Plus[i+1][j+1] * Minus[i+1][j+1] + Epsilon) / (pow (Plus[i+1][j+1], 2) + pow (Minus[i+1][j+1], 2) + Epsilon);
      EmptySetR[i][j] = EmptySet[i+1][j+1] - 0.25*s * ((1-kappa*s) * Plus[i+1][j+1] + (1+kappa*s) * Minus[i+1][j+1]); 
    }
  }


}




inline double CFD2D::RhoFluxXi (int i, int j) {
  return 0.5 * (RhoUR[i][j] + RhoUL[i][j] - RhoVR[i][j] - RhoVL[i][j]
		- (AbsJacobiA[0][0] * (RhoR[i][j] - RhoL[i][j]) 
		   + AbsJacobiA[0][1] * (RhoUR[i][j] - RhoUL[i][j]) 
		   + AbsJacobiA[0][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[0][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoUFluxXi (int i, int j) {
  return 0.5 * (RhoUR[i][j] * (UR[i][j] - VR[i][j]) + PR[i][j] 
		+ RhoUL[i][j] * (UL[i][j] - VL[i][j]) + PL[i][j]
		- (AbsJacobiA[1][0] * (RhoR[i][j] - RhoL[i][j]) 
		   + AbsJacobiA[1][1] * (RhoUR[i][j] - RhoUL[i][j]) 
		   + AbsJacobiA[1][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[1][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoVFluxXi (int i, int j) {
  return 0.5 * (RhoVR[i][j] * (UR[i][j] - VR[i][j]) - PL[i][j] 
		+ RhoVL[i][j] + (UL[i][j] - VL[i][j]) - PR[i][j]
		- (AbsJacobiA[2][0] * (RhoR[i][j] - RhoL[i][j])
		   + AbsJacobiA[2][1] * (RhoUR[i][j] - RhoUL[i][j])
		   + AbsJacobiA[2][2] * (RhoVR[i][j] - RhoVL[i][j])
		   + AbsJacobiA[2][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}
inline double CFD2D::RhoEFluxXi (int i, int j) {
  return 0.5 * ((RhoUR[i][j] - RhoVR[i][j]) * HR[i][j]
		+ (RhoUL[i][j] - RhoVL[i][j]) * HL[i][j]
		- (AbsJacobiA[3][0] * (RhoR[i][j] - RhoL[i][j])
		   + AbsJacobiA[3][1] * (RhoUR[i][j] - RhoUL[i][j])
		   + AbsJacobiA[3][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[3][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}

void CFD2D::AbsLambdaXi (int i, int j) {

  AbsLambda[0][0] = fabs (RoeU[i][j] - RoeV[i][j]) * RI2 - RoeC[i][j];   AbsLambda[0][1] = 0.0;   AbsLambda[0][2] = 0.0;  AbsLambda[0][3] = 0.0;
  AbsLambda[1][0] = 0.0;   AbsLambda[1][1] = fabs (RoeU[i][j] - RoeV[i][j]) * RI2;   AbsLambda[1][2] = 0.0;  AbsLambda[1][3] = 0.0;
  AbsLambda[2][0] = 0.0;   AbsLambda[2][2] = 0.0;   AbsLambda[2][2] = fabs (RoeU[i][j] - RoeV[i][j]) * RI2 + RoeC[i][j];  AbsLambda[2][3] = 0.0;
  AbsLambda[3][0] = 0.0;   AbsLambda[3][1] = 0.0;   AbsLambda[3][2] = 0.0;  AbsLambda[3][3] = fabs (RoeU[i][j] - RoeV[i][j]) * RI2;
}


inline void CFD2D::RightDiagonalMatrixXi (int i, int j) {
//     int ii, kk; 

  //    for (j = 0; j < 3; j++)
  //      DiagR[j][0] = 1.0;
  DiagR[0][0] = 1.0;
  DiagR[0][1] = 1.0;
  DiagR[0][2] = 1.0;
  DiagR[0][3] = 0.0;

  DiagR[1][0] = RoeU[i][j] - RI2 * RoeC[i][j];   
  DiagR[1][1] = RoeU[i][j];   
  DiagR[1][2] = RoeU[i][j] + RI2 * RoeC[i][j];  
  DiagR[1][3] = Rho[i][j];

  DiagR[2][0] = RoeV[i][j] + RI2 * RoeC[i][j];    
  DiagR[2][1] = RoeV[i][j];    
  DiagR[2][2] = RoeV[i][j] - RI2 * RoeC[i][j];    
  DiagR[2][3] = Rho[i][j];  

  DiagR[3][0] = RoeH[i][j] - RI2 * RoeC[i][j] * (RoeU[i][j] - RoeV[i][j]);
  DiagR[3][1] = 0.5 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));
  DiagR[3][2] = RoeH[i][j] + RI2 * RoeC[i][j] * (RoeU[i][j] - RoeV[i][j]);
  DiagR[3][3] = RoeRho[i][j] * (RoeU[i][j] + RoeV[i][j]);    

//    for (ii = 0; ii < 3; ii++) {
//      for (kk = 0; kk < 3; kk++)  printf ("%lf\t ", DiagR[ii][kk]);
//      printf ("\n");
//    }
//     printf ("\n");

}



inline void CFD2D::LeftDiagonalMatrixXi (int i, int j) {
  //  double RI2, RhoI2;

    //  Compose Left Diagonal Matrix for F|i-1/2
  //  Eta2 = Eta1 * (pow (RoeU, 2) + pow (RoeV, 2));

  RI2 = 1.0 / sqrt (2.0);
  RhoI2 = 0.5 * Rho[i][j];
  Eta2 = (Gamma - 1.0) / pow (RoeC[i][j], 2);
  Eta1 = 0.5 * Eta2 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));

  DiagL[0][0] = 0.5 * (Eta1 + RI2 * (RoeU[i][j] - RoeV[i][j])/ RoeC[i][j]);
  DiagL[0][1] = -0.5 * (Eta2 * RoeU[i][j] - RI2 / RoeC[i][j]);
  DiagL[0][2] = -0.5 * (Eta2 * RoeV[i][j] + RI2 / RoeC[i][j]);
  DiagL[0][3] = 0.5 * Eta2;

  DiagL[1][0] = 1.0 - Eta1;
  DiagL[1][1] = Eta2 * RoeU[i][j];
  DiagL[1][2] = Eta2 * RoeV[i][j];
  DiagL[1][3] = -Eta2;

  DiagL[2][0] = 0.5 * (Eta1 - (RoeU[i][j] - RoeV[i][j]) / RoeC[i][j]);
  DiagL[2][1] = -0.5 * (Eta2 * RoeU[i][j] - RI2 / RoeC[i][j]);
  DiagL[2][2] = -0.5 * (Eta2 * RoeV[i][j] + RI2 / RoeC[i][j]);
  DiagL[2][3] = 0.5 * Eta2;

  DiagL[3][0] = -(RoeU[i][j] + RoeV[i][j]) * RhoI2;
  DiagL[3][1] = RhoI2;
  DiagL[3][2] = RhoI2;
  DiagL[3][3] = 0.0;
}





//////////////////////////////////////////////////
//     Multiple Eta-Directional Differential    //
//////////////////////////////////////////////////
inline void CFD2D::EtaDirectionalFlux () {
  for (i = 2; i < LLX-2; i++) {
    for (j = 2; j < LLY-2; j++) {
      //      QLeftRightEta ();   //  Obtain RhoR, RhoL, RhoUR, RhoL, RhoER, RhoEL
      //
      //      RoeAverage ();

      //      Obtain Flux of X-directional Emptyset|i+1/2
      RightDiagonalMatrixEta (i, j);
      LeftDiagonalMatrixEta (i, j);
      //      printf ("Diagonal Matrix for X direction obtained\n");

      AbsLambdaEta (i, j);
      AbsJacobiAXY ();
      //      printf ("Jacobi matrix obtained X-direction obtained\n");

      RhoOutF = RhoFluxEta (i, j);
      RhoUOutF = RhoUFluxEta (i, j);
      RhoVOutF = RhoVFluxEta (i, j);
      RhoEOutF = RhoEFluxEta (i, j);

      //      Obtain Flux of X-directional Emptyset|i-1/2
      //      RoeAverage (i-1, j);
      //      RoeAverage ();
      RightDiagonalMatrixEta (i+1, j-1);
      LeftDiagonalMatrixEta (i+1, j-1);
      AbsLambdaEta (i+1, j-1);
      AbsJacobiAXY ();

      RhoInF = RhoFluxEta (i+1, j-1);
      RhoUInF = RhoUFluxEta (i+1, j-1);
      RhoVInF = RhoVFluxEta (i+1, j-1);
      RhoEInF = RhoEFluxEta (i+1, j-1);


      RhoFluxF[i][j]  = RhoOutF - RhoInF;
      RhoUFluxF[i][j] = RhoUOutF - RhoUInF;
      RhoVFluxF[i][j] = RhoVOutF - RhoVInF;
      RhoEFluxF[i][j] = RhoEOutF - RhoEInF;
    }
  }
}
inline void CFD2D::QLeftRightEta () {
   EmptySetEta (Rho, DeltaPlus, DeltaMinus, RhoL, RhoR);
   EmptySetEta (U, DeltaPlus, DeltaMinus, UL, UR);
   EmptySetEta (V, DeltaPlus, DeltaMinus, VL, VR);
   EmptySetEta (P, DeltaPlus, DeltaMinus, PL, PR);
   EmptySetEta (H, DeltaPlus, DeltaMinus, HL, HR);

   //  Rho|i+1/2, U|i+1/2, V|i+1/2 --> RhoUR, RhoUL, RhoVR, RhoVL
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoUL[i][j] = RhoL[i][j] * UL[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoUR[i][j] = RhoR[i][j] * UR[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoVL[i][j] = RhoL[i][j] * VL[i][j];
     }
   }
   for (i = 0; i < LLX; i++) {
     for (j = 0; j < LLY; j++) {
       RhoVR[i][j] = RhoR[i][j] * VR[i][j];
     }
   }

  //   Rho|i+1/2, P|i+1/2, U|i+1/2, V|i+1/2 -->
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoEL[i][j] = PL[i][j] / (Gamma - 1.0) + 0.5 * RhoL[i][j] * (pow (UL[i][j], 2) + pow (VL[i][j], 2));
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoER[i][j] = PR[i][j] / (Gamma - 1.0) + 0.5 * RhoR[i][j] * (pow (UR[i][j], 2) + pow (VR[i][j], 2));
    }
  }
}

inline void CFD2D::EmptySetEta (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR) {
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      Plus[i][j]  = minmod (EmptySet[i-1][j+1] - EmptySet[i][j], b*(EmptySet[i][j] - EmptySet[i+1][j-1]));
    }
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      Minus[i][j] = minmod (EmptySet[i][j] - EmptySet[i+1][j-1], b*(EmptySet[i-1][j+1] - EmptySet[i][j]));
      //      printf ("%d, %d, Rho = %le, DeltaMinus = %le\n", 
      //	      i, j, In[i][j], DeltaMinus[i][j]);
    }
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      s = (2.0 * Plus[i][j] * Minus[i][j] + Epsilon) / (pow (Plus[i][j], 2) + pow (Minus[i][j], 2) + Epsilon);
      EmptySetL[i][j] = EmptySet[i][j] + 0.25*s * ((1-kappa*s)* Plus[i][j] + (1+kappa*s) * Minus[i][j]); 
      //      printf ("%d, %d, Rho = %le, RhoL = %le\n", 
      //      	      i, j, EmptySet[i][j], EmptySetL[i][j]);
    }
  }  
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++) {
      s = (2.0 * Plus[i+1][j-1] * Minus[i+1][j-1] + Epsilon) / (pow (Plus[i+1][j-1], 2) + pow (Minus[i+1][j-1], 2) + Epsilon);
      EmptySetR[i][j] = EmptySet[i+1][j-1] - 0.25*s * ((1-kappa*s) * Plus[i+1][j-1] + (1+kappa*s) * Minus[i+1][j-1]); 
    }
  }


}




inline double CFD2D::RhoFluxEta (int i, int j) {
  return 0.5 * (RhoUR[i][j] + RhoUL[i][j] + RhoVR[i][j] + RhoVL[i][j]
		- (AbsJacobiA[0][0] * (RhoR[i][j] - RhoL[i][j]) 
		   + AbsJacobiA[0][1] * (RhoUR[i][j] - RhoUL[i][j]) 
		   + AbsJacobiA[0][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[0][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoUFluxEta (int i, int j) {
  return 0.5 * (RhoUR[i][j] * (UR[i][j] + VR[i][j]) + PR[i][j] 
		+ RhoUL[i][j] * (UL[i][j] + VL[i][j]) + PL[i][j]
		- (AbsJacobiA[1][0] * (RhoR[i][j] - RhoL[i][j]) 
		   + AbsJacobiA[1][1] * (RhoUR[i][j] - RhoUL[i][j]) 
		   + AbsJacobiA[1][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[1][3] * (RhoER[i][j] - RhoEL[i][j])));
}
inline double CFD2D::RhoVFluxEta (int i, int j) {
  return 0.5 * (RhoVR[i][j] * (VR[i][j] + UR[i][j]) + PL[i][j] 
		+ RhoVL[i][j] + (VL[i][j] + UL[i][j]) + PR[i][j]
		- (AbsJacobiA[2][0] * (RhoR[i][j] - RhoL[i][j])
		   + AbsJacobiA[2][1] * (RhoUR[i][j] - RhoUL[i][j])
		   + AbsJacobiA[2][2] * (RhoVR[i][j] - RhoVL[i][j])
		   + AbsJacobiA[2][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}
inline double CFD2D::RhoEFluxEta (int i, int j) {
  return 0.5 * ((RhoUR[i][j] + RhoVR[i][j]) * HR[i][j]
		+ (RhoUL[i][j] + RhoVL[i][j]) * HL[i][j]
		- (AbsJacobiA[3][0] * (RhoR[i][j] - RhoL[i][j])
		   + AbsJacobiA[3][1] * (RhoUR[i][j] - RhoUL[i][j])
		   + AbsJacobiA[3][2] * (RhoVR[i][j] - RhoVL[i][j]) 
		   + AbsJacobiA[3][3] * (RhoER[i][j] - RhoEL[i][j]))); 
}


void CFD2D::AbsLambdaEta (int i, int j) {

  AbsLambda[0][0] = fabs (RoeU[i][j] + RoeV[i][j]) * RI2 - RoeC[i][j];   AbsLambda[0][1] = 0.0;   AbsLambda[0][2] = 0.0;  AbsLambda[0][3] = 0.0;
  AbsLambda[1][0] = 0.0;   AbsLambda[1][1] = fabs (RoeU[i][j] + RoeV[i][j]) * RI2;   AbsLambda[1][2] = 0.0;  AbsLambda[1][3] = 0.0;
  AbsLambda[2][0] = 0.0;   AbsLambda[2][2] = 0.0;   AbsLambda[2][2] = fabs (RoeU[i][j] + RoeV[i][j]) * RI2 + RoeC[i][j];  AbsLambda[2][3] = 0.0;
  AbsLambda[3][0] = 0.0;   AbsLambda[3][1] = 0.0;   AbsLambda[3][2] = 0.0;  AbsLambda[3][3] = fabs (RoeU[i][j] + RoeV[i][j]) * RI2;
}


inline void CFD2D::RightDiagonalMatrixEta (int i, int j) {
//     int ii, kk; 

  //    for (j = 0; j < 3; j++)
  //      DiagR[j][0] = 1.0;
  DiagR[0][0] = 1.0;
  DiagR[0][1] = 1.0;
  DiagR[0][2] = 1.0;
  DiagR[0][3] = 0.0;

  DiagR[1][0] = RoeU[i][j] - RI2 * RoeC[i][j];   
  DiagR[1][1] = RoeU[i][j];   
  DiagR[1][2] = RoeU[i][j] + RI2 * RoeC[i][j];  
  DiagR[1][3] = RoeRho[i][j];

  DiagR[2][0] = RoeV[i][j] - RI2 * RoeC[i][j];    
  DiagR[2][1] = RoeV[i][j];    
  DiagR[2][2] = RoeV[i][j] + RI2 * RoeC[i][j];    
  DiagR[2][3] = -RoeRho[i][j];  

  DiagR[3][0] = RoeH[i][j] - RI2 * RoeC[i][j] * (RoeU[i][j] + RoeV[i][j]);
  DiagR[3][1] = 0.5 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));
  DiagR[3][2] = RoeH[i][j] + RI2 * RoeC[i][j] * (RoeU[i][j] + RoeV[i][j]);
  DiagR[3][3] = RoeRho[i][j] * (RoeU[i][j] - RoeV[i][j]);    

//    for (ii = 0; ii < 3; ii++) {
//      for (kk = 0; kk < 3; kk++)  printf ("%lf\t ", DiagR[ii][kk]);
//      printf ("\n");
//    }
//     printf ("\n");

}



inline void CFD2D::LeftDiagonalMatrixEta (int i, int j) {
  //  double RI2, RhoI2;
    //  Compose Left Diagonal Matrix for F|i-1/2
  //  Eta2 = Eta1 * (pow (RoeU, 2) + pow (RoeV, 2));

  RI2 = 1.0 / sqrt (2.0);
  RhoI2 = 0.5 * Rho[i][j];
  Eta2 = (Gamma - 1.0) / pow (RoeC[i][j], 2);
  Eta1 = 0.5 * Eta2 * (pow (RoeU[i][j], 2) + pow (RoeV[i][j], 2));

  DiagL[0][0] = 0.5 * (Eta1 + RI2 * (RoeU[i][j] + RoeV[i][j]) / RoeC[i][j]);
  DiagL[0][1] = -0.5 * (Eta2 * RoeU[i][j] + RI2 / RoeC[i][j]);
  DiagL[0][2] = -0.5 * (Eta2 * RoeV[i][j] + RI2 / RoeC[i][j]);
  DiagL[0][3] = 0.5 * Eta2;

  DiagL[1][0] = 1.0 - Eta1;
  DiagL[1][1] = Eta2 * RoeU[i][j];
  DiagL[1][2] = Eta2 * RoeV[i][j];
  DiagL[1][3] = -Eta2;

  DiagL[2][0] = 0.5 * (Eta1 - RI2 * (RoeU[i][j] + RoeV[i][j]) / RoeC[i][j]);
  DiagL[2][1] = -0.5 * (Eta2 * RoeU[i][j] - RI2 / RoeC[i][j]);
  DiagL[2][2] = -0.5 * (Eta2 * RoeV[i][j] - RI2 / RoeC[i][j]);
  DiagL[2][3] = 0.5 * Eta2;

  DiagL[3][0] = -0.5 * (RoeU[i][j] - RoeV[i][j]) / RoeRho[i][j];
  DiagL[3][1] = 0.5 / RoeRho[i][j];
  DiagL[3][2] = -0.5 / RoeRho[i][j];
  DiagL[3][3] = 0.0;
}







//////////////////////////////////////////////////
//        LES Method for Turbulent Flow         //
//////////////////////////////////////////////////
void CFD2D::LES (int n, double dt) {
  rx = 1.0 / dx;
  ry = 1.0 / dy;
  Cs = 0.173;
  dLES = sqrt (dx*dy);
  rx2 = 1.0 / pow (dx, 2);
  ry2 = 1.0 / pow (dy, 2);
  rxry = rx * ry;
  Theta = 0.0;
  Mu = 1.8e-05;

  //  printf ("hoge\n");

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)
      MuE[i][j] = Mu + 0.5 * Rho[i][j] * Cs * dLES *
	0.25 * (sqrt (2.0 * (pow (rx * (U[i+1][j] - U[i-1][j]), 2) + 
			     pow (ry * (V[i][j+1] - V[i][j-1]), 2) +
			     rx*ry* (U[i][j+1] - U[i][j-1]) * (V[i+1][j] - V[i-1][j])) +
		      pow (rx * (V[i+1][j] - V[i-1][j]), 2) + 
		      pow (ry * (U[i][j+1] - U[i][j-1]), 2) 
 		      ));
//       MuE[i][j] = Mu * 0.5 * Rho[i][j] * Cs * dLES *
// 	0.25 * (sqrt (2.0 * (pow (rx * (U[i+1][j] - U[i-1][j]), 2) + 
// 			     pow (ry * (V[i][j+1] - V[i][j-1]), 2) +
// 			     rx*ry* (U[i][j+1] - U[i][j-1]) * (V[i+1][j] - U[i-1][j])) +
// 		      pow (rx * (V[i+1][j] - V[i-1][j]), 2) + 
// 		      pow (ry * (U[i][j+1] - U[i][j-1]), 2) 
//  		      ));
  }


  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)   RhoS[i][j] = 0.0;
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)
      //      RhoUS[i][j] = 0.25 * (MuE[i][j] * (U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) * rx2 
      //			    + (Theta + 1.0 / 3.0 * MuE[i][j]) * 
      //			    ((U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) * rx2 + (V[i+1][j+1] - V[i-1][j-1]) * rxry)) * dt;
      RhoUS[i][j] = (MuE[i][j] * ((U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) * rx2
				  + (U[i][j+1] - 2.0 * U[i][j] + U[i][j+1]) * ry2)
		     + (Theta + 1.0 / 3.0 * MuE[i][j]) * 
		     ((U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) * rx2
		      + 0.25 * (V[i+1][j+1] - V[i-1][j-1]) * rxry)) * dt;
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)
      //      RhoVS[i][j] = 0.25 * (MuE[i][j] * (V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) * ry2 
      //			    + (Theta + 1.0 / 3.0 * MuE[i][j]) * 
      //			    ((V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) * ry2 + (U[i+1][j+1] - U[i-1][j-1]) * rxry)) * dt;
      RhoVS[i][j] = (MuE[i][j] * ((V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) * ry2 
				  + (V[i+1][j] - 2.0 * V[i][j] + V[i-1][j]) * rx2)
		     + (Theta + 1.0 / 3.0 * MuE[i][j]) * 
		     ((V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) * ry2 
		      + 0.25 * (U[i+1][j+1] - U[i-1][j-1]) * rxry)) * dt;
  }
  //  printf ("dt = %lf, rx = %lf, rxry = %lf\n", dt, rx, rxry);
  //  exit (1);

//   for (i = 1; i < LLX-1; i++) {
//     for (j = 1; i < LLY-1; j++)
//       SigmaXX[i][j] = 0.25 * (U[i+1][j] - 2.0 * U[i][j] + U[i-1][j]) * rx2;
//   }
//   for (i = 1; i < LLX-1; i++) {
//     for (j = 1; i < LLY-1; j++)
//       SmigaYY[i][j] = 0.25 * (V[i][j+1] - 2.0 * V[i][j] + V[i][j-1]) * ry2;
//   }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)
      TauXX[i][j] = 0.5 * (4.0/3.0 * (U[i+1][j] - U[i-1][j]) * rx - 2.0 / 3.0 * (V[i][j+1] - V[i][j-1]) * ry);
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)
      TauYY[i][j] = 0.5 * (4.0/3.0 * (V[i][j+1] - V[i][j-1]) * ry - 2.0 / 3.0 * (U[i+1][j] - U[i-1][j]) * rx);
  }
  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)
      TauXY[i][j] = 0.5 * ((U[i][j+1] - U[i][j-1]) * ry + (V[i+1][j] - V[i-1][j]) * rx);
  }

  for (i = 1; i < LLX-1; i++) {
    for (j = 1; j < LLY-1; j++)  
      RhoES[i][j] = 0.5 * MuE[i][j] * ((U[i+1][j] * TauXX[i+1][j] - U[i-1][j] * TauXX[i-1][j] 
				  + V[i+1][j] * TauXY[i+1][j] - V[i-1][j] * TauXY[i-1][j]) * rx
				 + (V[i][j+1] * TauYY[i][j+1] - V[i][j-1] * TauYY[i][j-1]
				    + U[i][j+1] * TauXY[i][j+1] - U[i][j-1] * TauXY[i][j-1]) * ry) * dt;
  }



}









//////////////////////////////////////////////////
//        Initial & Boundary Conditions         //
//////////////////////////////////////////////////
void CFD2D::InitFlat () {

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      Rho[i][j] = Rho0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      U[i][j] = U0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      V[i][j] = V0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
    P[i][j] = P0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
    RhoU[i][j] = Rho[i][j] * U[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoE[i][j] =  P[i][j] / (Gamma - 1.0) + 0.5 * Rho[i][j] * pow (U[i][j], 2);
    }
  }

}



void CFD2D::InitShock () {

  for (i = 0; i < LLX/2; i++) {
    for (j = 0; j < LLY; j++) {
      Rho[i][j] = 8.0 * Rho0;
//       if (Rho[i][j] < 1.0) {
// 	printf ("Rho inf %le\n", Rho[i][j]);
// 	exit (1);
//       }
    }
  }
  for (i = LLX/2; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      Rho[i][j] = Rho0;
    }
  }

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      U[i][j] = U0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      V[i][j] = V0;
    }
  }
  for (i = 0; i < LLX/2; i++) {
    for (j = 0; j < LLY; j++) {
    P[i][j] = 10.0 * P0;
    }
  }
  for (i = LLX/2; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
    P[i][j] = P0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
    RhoU[i][j] = Rho[i][j] * U[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoE[i][j] =  P[i][j] / (Gamma - 1.0) + 0.5 * Rho[i][j] * pow (U[i][j], 2);
    }
  }

}



void CFD2D::InitCaramel () {

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      Rho[i][j] = Rho0;
    }
  }


  for (i = 0; i < LLX/2; i++) {
    for (j = 0; j < LLY/2; j++) {
      Rho[i][j] = 8.0 * Rho0;
//       if (Rho[i][j] < 1.0) {
// 	printf ("Rho inf %le\n", Rho[i][j]);
// 	exit (1);
//       }
    }
  }

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      U[i][j] = U0;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      V[i][j] = V0;
    }
  }

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
    P[i][j] = P0;
    }
  }
  for (i = 0; i < LLX/2; i++) {
    for (j = 0; j < LLY/2; j++) {
    P[i][j] = 10.0 * P0;
    }
  }

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
    RhoU[i][j] = Rho[i][j] * U[i][j];
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      RhoE[i][j] =  P[i][j] / (Gamma - 1.0) + 0.5 * Rho[i][j] * pow (U[i][j], 2);
    }
  }

}



void CFD2D::Output1D (int j, int n, double t) {
  int i;

  sprintf (fuga, "out.1d.%03d", n);
  if (fp_out_1d = fopen (fuga, "w")) printf ("open %s for %dth output\n", fuga, n);
  else exit (1);
  printf ("outfreq = %d\n", output_freq);

  fprintf (fp_out, "#t = %lf\n", t);
  fprintf (fp_out, "#%d\t%d\t%lf\t%lf\n", LLX, LLY, Lx, Ly);

  for (i = 0; i < LLX; i++) {
    //    for (j = 0; j < LLY; j++) {
      //      fprintf (fp_out, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", 
      // 	       //         1    2    3    4    5    6    7    8   9    10
      fprintf (fp_out, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", 
       	       //         1    2    3    4    5    6    7    8    9  
      i*dx, Rho[i][j], U[i][j], V[i][j], P[i][j], RhoE[i][j], C[i][j], MuE[i][j],  sqrt(pow (U[i][j]/C[i][j], 2) + pow (V[i][j]/C[i][j], 2)));
    //  1        2       3          4       5       6           7       8             9        
      //    }
  }
  //  for (i = 0; i < LL; i++)    printf ("%d, %lf\n", i, Rho[i]);
  fclose (fp_out);
  //  printf ("hoge\n");
}


void CFD2D::Output (int n, double t) {
  int i;

  sprintf (fuga, "out.2d.%03d", n);
  if (fp_out = fopen (fuga, "w")) printf ("open %s for %dth output\n", fuga, n);
  else exit (1);
  printf ("outfreq = %d\n", output_freq);

  fprintf (fp_out, "#t = %lf\n", dt);
  fprintf (fp_out, "#%d\t%d\t%lf\t%lf\n", LLX, LLY, Lx, Ly);

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      //      fprintf (fp_out, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", 
      // 	       //         1    2    3    4    5    6    7    8   9    10
      fprintf (fp_out, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", 
       	       //         1    2    3    4    5    6    7    8    9   10
      i*dx, j*dy, Rho[i][j], U[i][j], V[i][j], P[i][j], RhoE[i][j], C[i][j], MuE[i][j],  sqrt(pow (U[i][j]/C[i][j], 2) + pow (V[i][j]/C[i][j], 2)));
    //  1    2        3          4       5       6           7       8             9        10
      // 	       i*dx, j*dy, Rho[i][j], U[i][j], V[i][j], P[i][j], RhoE[i][j], C[i][j],  RhoL[i][j], RhoR[i][j]);
      //         1    2        3          4         5      6           7            8             9        10
    //       	       i*dx, j*dy, Rho[i][j], U[i][j], P[i][j], RhoE[i][j], C[i][j],  RhoFluxE[i][j], RhoFluxF[i][j]);
      //         1    2        3          4         5      6           7            8             9
      //  for (i = 0; i < LL; i++)    printf ("%d, %lf\n", i, Rho[i]);
    }
  }
  //  for (i = 0; i < LL; i++)    printf ("%d, %lf\n", i, Rho[i]);
  fclose (fp_out);
  //  printf ("hoge\n");
}


//////////////////////////////////////////////////
//              Boundary Condition              //
//////////////////////////////////////////////////
void CFD2D::NeumannX () {
  for (j = 0; j < LLY; j++) {
    Rho[2][j] = Rho[3][j];   Rho[1][j] = Rho[2][j];   Rho[0][j] = Rho[1][j];
  }
  for (j = 0; j < LLY; j++) {
    Rho[LLX-3][j] = Rho[LLX-4][j];  Rho[LLX-2][j] = Rho[LLX-3][j];  Rho[LLX-1][j] = Rho[LLX-2][j];
  }
  for (j = 0; j < LLY; j++) {
    RhoU[2][j] = RhoU[3][j];   RhoU[1][j] = RhoU[2][j];   RhoU[0][j] = RhoU[1][j];
  }
  for (j = 0; j < LLY; j++) {
    RhoU[LLX-3][j] = RhoU[LLX-4][j];  RhoU[LLX-2][j] = RhoU[LLX-3][j];  RhoU[LLX-1][j] = RhoU[LLX-2][j];
  }
  for (j = 0; j < LLY; j++) {
    RhoV[2][j] = RhoV[3][j];   RhoV[1][j] = RhoV[2][j];   RhoV[0][j] = RhoV[1][j];
  }
  for (j = 0; j < LLY; j++) {
    RhoV[LLX-3][j] = RhoV[LLX-4][j];  RhoV[LLX-2][j] = RhoV[LLX-3][j];  RhoV[LLX-1][j] = RhoV[LLX-2][j];
  }
  for (j = 0; j < LLY; j++) {
    RhoE[2][j] = RhoE[3][j];   RhoE[1][j] = RhoE[2][j];   RhoE[0][j] = RhoE[1][j];
  }
  for (j = 0; j < LLY; j++) {
    RhoE[LLX-3][j] = RhoE[LLX-4][j];  RhoE[LLX-2][j] = RhoE[LLX-3][j];  RhoE[LLX-1][j] = RhoE[LLX-2][j];
  }

  for (j = 0; j < LLY; j++) {
    U[2][j] = U[3][j];   U[1][j] = U[2][j];   U[0][j] = U[1][j];
  }
  for (j = 0; j < LLY; j++) {
    U[LLX-3][j] = U[LLX-4][j];  U[LLX-2][j] = U[LLX-3][j];  U[LLX-1][j] = U[LLX-2][j];
  }


  for (j = 0; j < LLY; j++) {
    V[2][j] = V[3][j];   V[1][j] = V[2][j];   V[0][j] = V[1][j];
  }
  for (j = 0; j < LLY; j++) {
    V[LLX-3][j] = V[LLX-4][j];  V[LLX-2][j] = V[LLX-3][j];  V[LLX-1][j] = V[LLX-2][j];
  }

  for (j = 0; j < LLY; j++) {
    P[2][j] = P[3][j];   P[1][j] = P[2][j];   P[0][j] = P[1][j];
  }
  for (j = 0; j < LLY; j++) {
    P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
  }
}

void CFD2D::NeumannY () {
  for (i = 0; i < LLX; i++) {
      Rho[i][2] = Rho[i][3];      Rho[i][1] = Rho[i][2];      Rho[i][0] = Rho[i][1];
  }
  for (i = 0; i < LLX; i++) {
      Rho[i][LLY-3] = Rho[i][LLY-4];      Rho[i][LLY-2] = Rho[i][LLY-3];      Rho[i][LLY-1] = Rho[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
      RhoU[i][2] = RhoU[i][3];      RhoU[i][1] = RhoU[i][2];      RhoU[i][0] = RhoU[i][1];
  }
  for (i = 0; i < LLX; i++) {
      RhoU[i][LLY-3] = RhoU[i][LLY-4];      RhoU[i][LLY-2] = RhoU[i][LLY-3];      RhoU[i][LLY-1] = RhoU[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
      RhoV[i][2] = RhoV[i][3];      RhoV[i][1] = RhoV[i][2];      RhoV[i][0] = RhoV[i][1];
  }
  for (i = 0; i < LLX; i++) {
    RhoV[i][LLY-3] = RhoV[i][LLY-4];      RhoV[i][LLY-2] = RhoV[i][LLY-3];      RhoV[i][LLY-1] = RhoV[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
      RhoE[i][2] = RhoE[i][3];      RhoE[i][1] = RhoE[i][2];      RhoE[i][0] = RhoE[i][1];
  }
  for (i = 0; i < LLX; i++) {
      RhoE[i][LLY-3] = RhoE[i][LLY-4];      RhoE[i][LLY-2] = RhoE[i][LLY-3];      RhoE[i][LLY-1] = RhoE[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
      U[i][2] = U[i][3];      U[i][1] = U[i][2];      U[i][0] = U[i][1];
  }
  for (i = 0; i < LLX; i++) {
      U[i][LLY-3] = U[i][LLY-4];      U[i][LLY-2] = U[i][LLY-3];      U[i][LLY-1] = U[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
      V[i][2] = V[i][3];      V[i][1] = V[i][2];      V[i][0] = V[i][1];
  }
  for (i = 0; i < LLX; i++) {
      V[i][LLY-3] = V[i][LLY-4];      V[i][LLY-2] = V[i][LLY-3];      V[i][LLY-1] = V[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
      P[i][2] = P[i][3];      P[i][1] = P[i][2];      P[i][0] = P[i][1];
  }
  for (i = 0; i < LLX; i++) {
      P[i][LLY-3] = P[i][LLY-4];      P[i][LLY-2] = P[i][LLY-3];      P[i][LLY-1] = P[i][LLY-2];
  }
}



void CFD2D::NoSlipX () {
  for (j = 0; j < LLY; j++) {
    Rho[2][j] = Rho[3][j];   Rho[1][j] = Rho[2][j];   Rho[0][j] = Rho[1][j];
  }

  for (j = 0; j < LLY; j++) {
    Rho[LLX-3][j] = Rho[LLX-4][j];  Rho[LLX-2][j] = Rho[LLX-3][j];  Rho[LLX-1][j] = Rho[LLX-2][j];
  }

  for (j = 0; j < LLY; j++) {
    U[2][j] = 0.0;   U[1][j] = 0.0;   U[0][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    U[LLX-3][j] = 0.0;  U[LLX-2][j] = 0.0;  U[LLX-1][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    V[2][j] = 0.0;   V[1][j] = 0.0;   V[0][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    V[LLX-3][j] = 0.0;  V[LLX-2][j] = 0.0;  V[LLX-1][j] = 0.0;
  }



  for (j = 0; j < LLY; j++) {
    RhoU[2][j] = 0.0;   RhoU[1][j] = 0.0;   RhoU[0][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    RhoU[LLX-3][j] = 0.0;  RhoU[LLX-2][j] = 0.0;  RhoU[LLX-1][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    RhoV[2][j] = 0.0;   RhoV[1][j] = 0.0;   RhoV[0][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    RhoV[LLX-3][j] = 0.0;  RhoV[LLX-2][j] = 0.0;  RhoV[LLX-1][j] = 0.0;
  }


  for (j = 0; j < LLY; j++) {
    P[2][j] = P[3][j];   P[1][j] = P[2][j];   P[0][j] = P[1][j];
  }
  for (j = 0; j < LLY; j++) {
    P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
  }

  for (j = 0; j < LLY; j++) {
    RhoE[2][j] = P[2][j] / (Gamma - 1.0);
    RhoE[1][j] = P[1][j] / (Gamma - 1.0);
    RhoE[0][j] = P[0][j] / (Gamma - 1.0);
  }
  for (j = 0; j < LLY; j++) {
    RhoE[LLX-3][j] = P[LLX-3][j] / (Gamma - 1.0);
    RhoE[LLX-2][j] = P[LLX-2][j] / (Gamma - 1.0);
    RhoE[LLX-1][j] = P[LLX-1][j] / (Gamma - 1.0);
  }


}




void CFD2D::NoSlipY () {

  for (i = 0; i < LLX; i++) {
      Rho[i][2] = Rho[i][3];      Rho[i][1] = Rho[i][2];      Rho[i][0] = Rho[i][1];
  }
  for (i = 0; i < LLX; i++) {
      Rho[i][LLY-3] = Rho[i][LLY-4];      Rho[i][LLY-2] = Rho[i][LLY-3];      Rho[i][LLY-1] = Rho[i][LLY-2];
  }


  for (i = 0; i < LLX; i++) {
      U[i][2] = 0.0;      U[i][1] = 0.0;      U[i][0] = 0.0;
  }
  for (i = 0; i < LLX; i++) {
      U[i][LLY-3] = 0.0;      U[i][LLY-2] = 0.0;      U[i][LLY-1] = 0.0;
  }

  for (i = 0; i < LLX; i++) {
      V[i][2] = 0.0;      V[i][1] = 0.0;      V[i][0] = 0.0;
  }
  for (i = 0; i < LLX; i++) {
      V[i][LLY-3] = 0.0;      V[i][LLY-2] = 0.0;      V[i][LLY-1] = 0.0;
  }

  for (i = 0; i < LLX; i++) {
      RhoU[i][2] = 0.0;      RhoU[i][1] = 0.0;      RhoU[i][0] = 0.0;
  }
  for (i = 0; i < LLX; i++) {
      RhoU[i][LLY-3] = 0.0;      RhoU[i][LLY-2] = 0.0;      RhoU[i][LLY-1] = 0.0;
  }

  for (i = 0; i < LLX; i++) {
      RhoV[i][2] = 0.0;      RhoV[i][1] = 0.0;      RhoV[i][0] = 0.0;
  }
  for (i = 0; i < LLX; i++) {
    RhoV[i][LLY-3] = 0.0;      RhoV[i][LLY-2] = 0.0;      RhoV[i][LLY-1] = 0.0;
  }


  for (i = 0; i < LLX; i++) {
      P[i][2] = P[i][3];      P[i][1] = P[i][2];      P[i][0] = P[i][1];
  }
  for (i = 0; i < LLX; i++) {
      P[i][LLY-3] = P[i][LLY-4];      P[i][LLY-2] = P[i][LLY-3];      P[i][LLY-1] = P[i][LLY-2];
  }

  for (i = 0; i < LLX; i++) {
    RhoE[i][2] = P[i][2] / (Gamma - 1.0);
    RhoE[i][1] = P[i][1] / (Gamma - 1.0);
    RhoE[i][0] = P[i][0] / (Gamma - 1.0);
  }
  for (i = 0; i < LLX; i++) {
    RhoE[i][LLY-3] = P[i][LLY-3] / (Gamma - 1.0);
    RhoE[i][LLY-2] = P[i][LLY-2] / (Gamma - 1.0);
    RhoE[i][LLY-1] = P[i][LLY-1] / (Gamma - 1.0);
  }


}


void CFD2D::InletX () {
  for (j = 0; j < LLY; j++) {
    //    Rho[2][j] = Rho[3][j];   Rho[1][j] = Rho[2][j];   Rho[0][j] = Rho[1][j];
    Rho[2][j] = Rho0;   Rho[1][j] = Rho0;   Rho[0][j] = Rho0;
  }
  for (j = 0; j < LLY; j++) {
    U[2][j] = U[3][j];   U[1][j] = U[2][j];   U[0][j] = U[1][j];
  }
  for (j = 0; j < LLY; j++) {
    V[2][j] = V[3][j];   V[1][j] = V[2][j];   V[0][j] = V[1][j];
  }
  for (j = 0; j < LLY; j++) {
    P[2][j] = P0;   P[1][j] = P0;   P[0][j] = P0;
  }

}


void CFD2D::OutletX (double t, double dt) {

  for (j = 0; j < LLY; j++) {
    if ((0.001 < t && t < 0.05) || (0.15 < t && t < 0.23)) {
      if (15 <= j && j <= 18) {
	Rho[LLX-3][j] = 0.9 * Rho0;  Rho[LLX-2][j] = 0.9 * Rho0;  Rho[LLX-1][j] = 0.9 * Rho0;
	U[LLX-3][j] = U[LLX-4][j];  U[LLX-2][j] = U[LLX-3][j];  U[LLX-1][j] = U[LLX-2][j];
	V[LLX-3][j] = V[LLX-4][j];  V[LLX-2][j] = V[LLX-3][j];  V[LLX-1][j] = V[LLX-2][j];
	P[LLX-3][j] = 0.9 * P0;  P[LLX-2][j] = 0.9 * P0;  P[LLX-1][j] = 0.9 * P0;
      } else {
	Rho[LLX-3][j] = Rho[LLX-4][j];  Rho[LLX-2][j] = Rho[LLX-3][j];  Rho[LLX-1][j] = Rho[LLX-2][j];
	U[LLX-3][j] = 0.0;  U[LLX-2][j] = 0.0;  U[LLX-1][j] = 0.0;
	V[LLX-3][j] = 0.0;  V[LLX-2][j] = 0.0;  V[LLX-1][j] = 0.0;
	P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
      }
    } else {
      Rho[LLX-3][j] = Rho[LLX-4][j];  Rho[LLX-2][j] = Rho[LLX-3][j];  Rho[LLX-1][j] = Rho[LLX-2][j];
      U[LLX-3][j] = 0.0;  U[LLX-2][j] = 0.0;  U[LLX-1][j] = 0.0;
      V[LLX-3][j] = 0.0;  V[LLX-2][j] = 0.0;  V[LLX-1][j] = 0.0;
      P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
    }
  }


  // for (j = 0; j < LLY; j++) {
//     if (0.001 < t && t < 0.002) {
//       if (15 <= j && j <= 18) {
// 	U[LLX-3][j] = U[LLX-4][j];  U[LLX-2][j] = U[LLX-3][j];  U[LLX-1][j] = U[LLX-2][j];
//       } else {
// 	U[LLX-3][j] = 0.0;  U[LLX-2][j] = 0.0;  U[LLX-1][j] = 0.0;
//       }
//     } else {
//       U[LLX-3][j] = 0.0;  U[LLX-2][j] = 0.0;  U[LLX-1][j] = 0.0;
//     }
//   }


//   for (j = 0; j < LLY; j++) {
//     if (0.001 <= t && t <= 0.002) {
//       if (15 <= j && j <= 18) {
// 	V[LLX-3][j] = V[LLX-4][j];  V[LLX-2][j] = V[LLX-3][j];  V[LLX-1][j] = V[LLX-2][j];
//       } else {
// 	V[LLX-3][j] = 0.0;  V[LLX-2][j] = 0.0;  V[LLX-1][j] = 0.0;
//       }
//     } else {
//       V[LLX-3][j] = 0.0;  V[LLX-2][j] = 0.0;  V[LLX-1][j] = 0.0;
//     }
//   }


//   for (j = 0; j < LLY; j++) {
//     if (0.001 < t && t < 0.002) {
//       if (15 <= j && j <= 18) {
// 	P[LLX-3][j] = 0.9 * P0;  P[LLX-2][j] = 0.9 * P0;  P[LLX-1][j] = 0.9 * P0;

//       } else {
// 	P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
//       }
//     } else {
//       P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
//     }
//  }
}



void CFD2D::StopBoundaryX () {
  for (j = 0; j < LLY; j++) {
    Rho[LLX-3][j] = Rho[LLX-4][j];  Rho[LLX-2][j] = Rho[LLX-3][j];  Rho[LLX-1][j] = Rho[LLX-2][j];
  }
  for (j = 0; j < LLY; j++) {
    U[LLX-3][j] = 0.0;  U[LLX-2][j] = 0.0;  U[LLX-1][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    V[LLX-3][j] = 0.0;  V[LLX-2][j] = 0.0;  V[LLX-1][j] = 0.0;
  }
  for (j = 0; j < LLY; j++) {
    P[LLX-3][j] = P[LLX-4][j];  P[LLX-2][j] = P[LLX-3][j];  P[LLX-1][j] = P[LLX-2][j];
  }



}





void CFD2D::FuelOutY (double t, double dt) {

  for (i = 0; i < LLX; i++) {
      Rho[i][2] = Rho[i][3];      Rho[i][1] = Rho[i][2];      Rho[i][0] = Rho[i][1];
  }
  for (i = 0; i < LLX; i++) {
      U[i][2] = 0.0;      U[i][1] = 0.0;      U[i][0] = 0.0;
  }
  for (i = 0; i < LLX; i++) {
      V[i][2] = 0.0;      V[i][1] = 0.0;      V[i][0] = 0.0;
  }
  for (i = 0; i < LLX; i++) {
      P[i][2] = P[i][3];      P[i][1] = P[i][2];      P[i][0] = P[i][1];
  }


  for (i = 0; i < LLX; i++) {
    if (0.14 * n_pulse + 0.001 < t && t < 0.14 * n_pulse + 0.047) {
      if ((i < LLX - (n_pulse % 6 * 60 + 5)) && ((LLX - (n_pulse % 6 * 60 + 8)) < i)) {
	//	printf ("i %d, t %lf\n", i, t);
	//	exit (1);
      	Rho[i][LLY-3] = 0.9 * Rho0;            	Rho[i][LLY-2] = 0.9 * Rho0;            	Rho[i][LLY-1] = 0.9 * Rho0;      
	U[i][LLY-3] = U[i][LLY-4];      U[i][LLY-2] = U[i][LLY-3];      U[i][LLY-1] = U[i][LLY-2];
	V[i][LLY-3] = V[i][LLY-4];      V[i][LLY-2] = V[i][LLY-3];      V[i][LLY-1] = V[i][LLY-2];
      	P[i][LLY-3] = 0.9 * P0;            	P[i][LLY-2] = 0.9 * P0;            	P[i][LLY-1] = 0.9 * P0;      
      } else {
	Rho[i][LLY-3] = Rho[i][LLY-4];      Rho[i][LLY-2] = Rho[i][LLY-3];      Rho[i][LLY-1] = Rho[i][LLY-2];
	U[i][LLY-3] = 0.0;      U[i][LLY-2] = 0.0;      U[i][LLY-1] = 0.0;
	V[i][LLY-3] = 0.0;      V[i][LLY-2] = 0.0;      V[i][LLY-1] = 0.0;
	P[i][LLY-3] = P[i][LLY-4];      P[i][LLY-2] = P[i][LLY-3];      P[i][LLY-1] = P[i][LLY-2];
      }
    } else {
      Rho[i][LLY-3] = Rho[i][LLY-4];      Rho[i][LLY-2] = Rho[i][LLY-3];      Rho[i][LLY-1] = Rho[i][LLY-2];
      U[i][LLY-3] = 0.0;      U[i][LLY-2] = 0.0;      U[i][LLY-1] = 0.0;
      V[i][LLY-3] = 0.0;      V[i][LLY-2] = 0.0;      V[i][LLY-1] = 0.0;
      P[i][LLY-3] = P[i][LLY-4];      P[i][LLY-2] = P[i][LLY-3];      P[i][LLY-1] = P[i][LLY-2];
    }
  }

  //    if (0.14 * n_pulse + 0.001 < t && t < 0.14 * n_pulse + 0.047) {
  if ((0.14 * n_pulse + 0.047 - t) * (0.14 * n_pulse + 0.047 - t - dt) < 0.0) {
    n_pulse++;
    printf ("pulse %d, t %lf\n", n_pulse, t);
    //    getch ();
    //    exit (1);
  }

}




//////////////////////////////////////////////////
//               Immersed Boundary              //
//////////////////////////////////////////////////
void CFD2D::LoadMask () {
  printf ("test to load mask\n");
  if (fp_mask = fopen ("inp.mask", "r")) {
    printf ("hoge\n");
  }

  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      fscanf (fp_mask, "%d\t%d\t%lf\n", &i_dummy, &j_dummy, &Mask[i][j]);
    }
  }

  fclose (fp_mask);
}

//void CFD2D::MaskMotion () {
//
//}


void CFD2D::MaskBoundaryCondition () {
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      if (0.9 < Mask[i][j]) Rho[i][j] = Rho0;
    }
  }       
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      if (0.9 < Mask[i][j]) U[i][j] = 0.0;
    }
  }         for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      if (0.9 < Mask[i][j]) V[i][j] = 0.0;
    }
  }       
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      if (0.9 < Mask[i][j]) P[i][j] = P0;
    }
  }       
}


void CFD2D::ImmersedBoundaryTrace () {
  int i_gc, j_gc, num_gc_mark;
  double AbsGradMask, rx, ry;


  rx = 0.5 / dx;
  ry = 0.5 / dy;

  fp_out = fopen ("out.mask.surface", "w");

  //  Searching for Points in the surface of Materials and Fluids
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      num_gc_mark = 0;
      for (i_gc = -1; i_gc < 2; i_gc++) {
	for (j_gc = -1; j_gc < 2; j_gc++) {
	  if ((Mask[i][j] - 0.5) * (Mask[i + i_gc][j + j_gc] - 0.5) < 0.0 && 0 < Mask[i][j]) {
	    num_gc_mark++;
	  }
	}
      }
      if (0 < num_gc_mark) {
	MaskSurface[i][j] = 1.0;
	//	fprintf (fp_out, "%d\t%d\t%d\n", i, j, num_gc_mark);
	//	fprintf (fp_out, "%d\t%d\t%d\n", i, j, 1);
      } else {
	MaskSurface[i][j] = 0.0;
	//	fprintf (fp_out, "%d\t%d\t%d\n", i, j, 0);
      }
    }
  }


  //  Searching for Ghost-Cell
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      MaskNormalX [i][j] = (Mask[i+1][j] - Mask[i-1][j]) / rx;
    }
  }
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      MaskNormalY [i][j] = (Mask[i][j+1] - Mask[i][j-1]) / rx;
    }
  }

  //  for (i = 0; i < LLX; i++) {
  //    for (j = 0; j < LLY; j++) {
      //      AbsGradMask = sqrt (pow (MaskNormalX[i][j], 2) + pow (MaskNormalY[i][j], 2));
      //      fprintf (fp_out, "%d\t%d\t%le\t%le\t%le\t%le\n", 
      //	       //        1   2    3     4   5    6
      //	       i, j, MaskSurface[i][j], MaskNormalX[i][j], MaskNormalY[i][j], AbsGradMask);
      //       1  2            3               4                     5             6
      //	       i, j, MaskSurface[i][j], MaskNormalX[i][j] / AbsGradMask, MaskNormalY[i][j] / AbsGradMask, AbsGradMask);
      //      //       1  2            3                      4                                5                      6
  //    }
  //  }


  //  fclose (fp_out);


}


void CFD2D::ImmersedBoundaryGhostCell () {
  int i_gc, j_gc, num_gc_mark, 
    near_ghost_cell_x[5], near_ghost_cell_y[5], num_nearby_cell;
  double NearByCellR[5];

  //  Searching for Points in the surface of Materials and Fluids
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      num_gc_mark = 0;
      for (i_gc = -1; i_gc < 2; i_gc++) {
	for (j_gc = -1; j_gc < 2; j_gc++) {
	  if ((Mask[i][j] - 0.5) * (Mask[i + i_gc][j + j_gc] - 0.5) < 0.0 && 0 < Mask[i][j]) {
	    num_gc_mark++;
	  }
	}
      }
      if (0 < num_gc_mark) {
	MaskSurface[i][j] = 1.0;
      } else {
	MaskSurface[i][j] = 0.0;
      }
    }
  }

  //  Searching for Nearby Cell to GhostCell
  for (i = 0; i < LLX; i++) {
    for (j = 0; j < LLY; j++) {
      num_nearby_cell = 0;
      if (0.5 < MaskSurface[i][j]) {
	for (im = i - 5; im < i + 5; im++) {
	  for (jm = j - 5; jm < j + 5; jm++) {
	    if (Mask[i][j] < 0.5) {
	      near_ghost_cell_x[num_nearby_cell] = im; 
	      near_ghost_cell_y[num_nearby_cell] = jm;
	      num_nearby_cell++;
	      
	    }
	  }
	}
	//  Sort
      }
    }
  }


}


//void CFD2D::ImmersedBoundaryNonSlip () {
//
//}
