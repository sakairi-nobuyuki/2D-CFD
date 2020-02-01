#include "functions.h"
#include <cstdio>


class CFD2D {
 protected:
  int i, j, k, l, im, jm, km, lm, i_dummy, j_dummy;
  char piyo[64], fuga[64], mud[1024];
  FILE *fp_in, *fp_out, *fp_out_1d, *fp_mask, *fp_mask_test;
  //  double **Rho, **RhoU, **RhoE, **RhoF, **RhoUF, **RhoEF,
  //    **U, **H, **C, **P, Rho0, U0, P0,
  double **RhoU, **RhoV, **RhoE, 
    **RhoFluxE, **RhoUFluxE, **RhoVFluxE, **RhoEFluxE,
    **RhoFluxF, **RhoUFluxF, **RhoVFluxF, **RhoEFluxF,
    **RoeRho, **RoeU, **RoeV, **RoeH, **RoeC,
    **H, **C, Rho0, U0, V0, P0,
    **Rhotmp, **RhoUtmp, **RhoVtmp, **RhoEtmp,
    **RhoR, **RhoL, **RhoUR, **RhoUL, **RhoVR, **RhoVL, **RhoEL, **RhoER, 
    **UR, **UL, **VR, **VL, **HR, **HL, **PR, **PL, **SlopeTmp, 
    //    RoeRho, RoeU, RoeV, RoeH, RoeC, Eta1, Eta2, q2, 
    Eta1, Eta2, q2, 
    RhoOutE, RhoUOutE, RhoVOutE, RhoEOutE, RhoInE, RhoUInE, RhoVInE, RhoEInE, 
    RhoOutF, RhoUOutF, RhoVOutF, RhoEOutF, RhoInF, RhoUInF, RhoVInF, RhoEInF, 
    RhoURtmp, RhoULtmp, RhoVRtmp, RhoVLtmp, RhoERtmp, RhoELtmp, HRtmp, HLtmp, 
    DiagR[4][4], DiagL[4][4], AbsLambda[4][4], AbsJacobiA[4][4], AbsJacobiAtmp[4][4],
    JacobiA[4][4], JacobiAtmp[4][4], Lambda[4][4], 
    RhoBtmp[3], UBtmp[3], PBtmp[3], RhoB[3], UB[3], PB[3], RhoUB[3], RhoUBtmp[3], RhoVB[3], RhoVBtmp[3], RhoEBtmp[3], RhoEB[3], 
    //  var. for 1 var. convection eq.
    **Q, **Qtmp, **E, **F, **DeltaPlus, **DeltaMinus, **QR, **QL, **FluxPlus, **FluxMinus,
    Lx, Ly, Co, s, kappa, b, Epsilon, dt, dx, dy, cx, cy, rx, ry, MaxPhaseSpeed, MaxU, MaxV, MaxC, MinP,
    // var. for LES
    **MuE, Mu, Theta, Cs, **RhoS, **RhoUS, **RhoVS, **RhoES, **TauXX, **TauYY, **TauXY, **SigmaXX, **SigmaYY,
    dLES, rx2, ry2, rxry,
    //  multi-direction derivertive
    **OblRhoFluxE, **OblRhoFluxF, RateOblique, RateNormal, dxi, deta, cxi, ceta, RI2, RhoI2,
    //  Immersed Boundary
    **Mask, **RhoFB, **RhoUFB, **RhoVFB, **RhoEFB, RefGC[5], **MaskSurface, **MaskNormalX, **MaskNormalY;
  //    NearGhostCellX[5], NearGhostCellY[5];
  inline void QLeftRightX ();
  inline void QLeftRightY ();
  //  inline void RoeAverage (int i, int j);
  inline void RoeAverage ();
  inline void RightDiagonalMatrix (int i, int j);
  inline void XDirectionalFlux ();
  inline void RightDiagonalMatrixX (int i, int j);
  inline void LeftDiagonalMatrixX (int i, int j);
  inline void AbsLambdaX (int i, int j);
  inline double RhoFluxX (int i, int j);
  inline double RhoUFluxX (int i, int j);
  inline double RhoVFluxX (int i, int j);
  inline double RhoEFluxX (int i, int j);
  inline void YDirectionalFlux ();
  inline void RightDiagonalMatrixY (int i, int j);
  inline void LeftDiagonalMatrixY (int i, int j);
  inline void AbsLambdaY (int i, int j);
  inline double RhoFluxY (int i, int j);
  inline double RhoUFluxY (int i, int j);
  inline double RhoVFluxY (int i, int j);
  inline double RhoEFluxY (int i, int j);
  inline void AbsJacobiAXY ();
  inline void DeltaPlusMnusX (double** In, double** DeltaPlus, double** DeltaMinus, double b);
  inline void EmptySetX (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR);
  inline void DeltaPlusMnusY (double** In, double** DeltaPlus, double** DeltaMinus, double b);
  inline void EmptySetLRY (double ** DeltaPlus, double** DeltaMinus, double** EmptySet, double** EmptySetL, double** EmptySetR);
  inline void EmptySetY (double **EmptySet, double ** Plus, double** Minus, double** EmptySetL, double** EmptySetR);
  inline double Slope (double DeltaPlus, double DeltaMinus, double Epsilon);
  //   Multiple Direction
  inline void XiDirectionalFlux ();
  inline void QLeftRightXi ();
  inline void EmptySetXi (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR);
  inline void DeltaPlusMnusXi (double** In, double** DeltaPlus, double** DeltaMinus, double b);
  inline void RightDiagonalMatrixXi (int i, int j);
  inline void LeftDiagonalMatrixXi (int i, int j);
  inline void AbsLambdaXi (int i, int j);
  inline double RhoFluxXi (int i, int j);
  inline double RhoUFluxXi (int i, int j);
  inline double RhoVFluxXi (int i, int j);
  inline double RhoEFluxXi (int i, int j);

  inline void EtaDirectionalFlux ();
  inline void QLeftRightEta ();
  //  inline void EmptySeEta (double** EmptySet, double** Plus, double** Minus, double** EmptySetL, double** EmptySetR);
  inline void EmptySetEta (double **EmptySet, double ** Plus, double** Minus, double** EmptySetL, double** EmptySetR);
  inline void DeltaPlusMnusEta (double** In, double** DeltaPlus, double** DeltaMinus, double b);
  inline void RightDiagonalMatrixEta (int i, int j);
  inline void LeftDiagonalMatrixEta (int i, int j);
  inline void AbsLambdaEta (int i, int j);
  inline double RhoFluxEta (int i, int j);
  inline double RhoUFluxEta (int i, int j);
  inline double RhoVFluxEta (int i, int j);
  inline double RhoEFluxEta (int i, int j);
  inline void CFL ();

 public:
  CFD2D ();
  CFD2D (char hoge[64]);
  ~CFD2D ();

  int output_freq, output_freq_1d, LLX, LLY, n_pulse;
  double Gamma, **Rho, **U, **V, **P;
  
  double RoeMUSCL (int n, double dt);
  void LES (int n, double dt);
  void InitShock ();
  void InitFlat ();
  void InitCaramel ();

  void NeumannX ();
  void NeumannY ();
  void NoSlipX ();
  void NoSlipY ();
  void FuelOutY (double t, double dt);
  void InletX ();
  void OutletX (double t, double dt);
  void StopBoundaryX ();

  void MaskBoundaryCondition ();
  void ImmersedBoundaryTrace ();
  void ImmersedBoundaryGhostCell ();
  //  void ImmersedBoundaryNonSlip ();
  void LoadMask ();
  void InterfaceToVolume (double RhoIn, double UIn, double PIn, double RhoOut, double UOut, double POut, int n);
  void xInterfaceToVolume (double RhoIn, double UIn, double PIn, double RhoOut, double UOut, double POut, int n);
  void xxInterfaceToVolume (double RhoIn, double UIn, double PIn, double RhoOut, double UOut, double POut, int n);
  void Output (int n, double t);
  void Output1D (int j, int n, double t);

};
