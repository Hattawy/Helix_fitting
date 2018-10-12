#include "math.h"
#include "stdio.h"
float eg6rtpc_eloss(int, float, float);

// Specs for EG6's Ne:DME (4:1) drift gas:
// DENSITY=1.03E-3; // density (g/cm^3)
// AEFF=126.787;    // effective atomic number (g/mol)
// ZEFF=66;         // effective atomic charge
// IEFF=99.794E-6;  // effective ionization energy (MeV)
// KAPPA=0.307075;  // PDG constant (MeV/g*cm^2)
// COEFF=KAPPA*ZEFF/AEFF*DENSITY // MeV/cm
// Scale factor to keep gain scaling near 1:
// COEFF100=COEFF*100
// COEFF800=COEFF*800
// COEFF1333=COEFF*1333
//#define COEFF100 0.01646459
//#define COEFF800 0.13171672
#define COEFF1333 0.219473 
#define IEFF  99.794E-6
#define MELEC 0.511
#define NPARTS 5

float eg6rtpc_dedx(float mom,float mass,float charge)
{ 
  // Bethe-Bloch from PDG.
  // "mom" is the real momentum in the drift region (NOT p/q).
  // "mom" and "mass" must be in units of MeV.
  // "charge" is unit charge.
  double beta,gamma,Tmax,logar,bb;
  beta  = mom/sqrt(mass*mass+mom*mom);
  gamma = 1./sqrt(1.-beta*beta);
  Tmax  = 2.*MELEC*pow(beta*gamma,2) / (1.+2.*gamma*MELEC/mass+pow(MELEC/mass,2));
  logar = log(2.*MELEC*pow(beta*gamma,2)*Tmax/pow(IEFF,2));
  bb    = COEFF1333 * pow(charge/beta,2) * (logar/2 - beta*beta);
  return (bb>0 ? bb : 0);
}

void eg6rtpc_pids(float praw,float dqdx,float theta,int* pid,float* preal)
{
  // input:
  // praw  = raw momentum in drift region (p/q) in units MeV
  // dqdx  = measured dE/dX in units ADC/mm
  // theta = polar angle in radians
  // output:
  // pid   = array of geant3 PIDs sorted by increasing distance to Bethe-Bloch
  // preal = array of energy-loss corrected vertex momentum corresponding to pid

  // The 5 hypotheses' (p,d,T,3He,4He) masses, charges, pids:
  const float mm[NPARTS]={938.272,1875.613,2808.921,2808.391,3727.379};
  const float qq[NPARTS]={1.0,1.0,1.0,2.0,2.0};
  const int geant3pids[NPARTS]={2212,45,46,49,47};
 
  int ii,itmp,swapped;
  float tmp,pdrift,ddqdx[NPARTS];

  for (ii=0; ii<NPARTS; ii++)
  {
    // unsorted pid:
    pid[ii] = ii;
    // Bethe-Bloch needs real momentum (NOT p/q):
    pdrift = praw*qq[ii];
    // distance to Bethe-Bloch for each hypothesis:
    ddqdx[ii]= fabs( dqdx - eg6rtpc_dedx(pdrift,mm[ii],qq[ii]) );
  }

  // bubble sort pid by ddqdx:
  swapped=1;
  while (swapped)
  {
      swapped=0;
      for (ii=1; ii<NPARTS; ii++)
      {
          if (ddqdx[ii] < ddqdx[ii-1])
          {
              tmp=ddqdx[ii];
              ddqdx[ii]=ddqdx[ii-1];
              ddqdx[ii-1]=tmp;

              itmp=pid[ii];
              pid[ii]=pid[ii-1];
              pid[ii-1]=itmp;
          
              swapped=1;
          }
      }
  }
 
  // apply energy loss corrections for the 5 hypotheses:
  for (ii=0; ii<NPARTS; ii++) 
  {
      // energy loss correction needs real momentum in drift region (NOT p/q):
      pdrift = praw*qq[pid[ii]];
      // convert to geant3 pid scheme (energy loss correction needs this too):
      pid[ii] = geant3pids[pid[ii]];
      // do energy loss correction:
      preal[ii] = eg6rtpc_eloss( pid[ii], theta, pdrift );
  }

}

