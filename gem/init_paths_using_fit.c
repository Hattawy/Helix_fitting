#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rtpc.h"

// last modified: M. Hattawy: 15/01/2013
// new parameter of delta_phi function
// fixing the parameters of the drift speed as a function of z 
//


// nab modified 2012
// - various problems fixed in combination with gem/clean_hits.c
// - hack for corner case
//

// These routines define the EG6 drift path calibration.
// Derived by fitting EG6 data.


void init_paths_using_fit_eg6()
{
  int ii, tt, num_bad=0, num_neg=0;
  float z;

  fprintf(stderr,"gem/init_paths_using_fit:  using fitted EG6 drift paths from pass 1 data\n");

  for (ii = 0; ii<NUM_PADS; ii++)
  {
    z = padLocs[ii].z;
    for (tt=0; tt<NAL_SAMP; tt++)
    {
      rawCYL[tt][ii].z = padLocs[ii].z;

      if(tt < TPC_TZERO) 
      {
        //rawCYL[tt][ii].s =SMIN; // technically this one is wrong
        //rawCYL[tt][ii].s =SMAX; // doesn't matter, hit will be rejected wether it's SMIN or SMAX
        rawCYL[tt][ii].s   = 9999.;
        rawCYL[tt][ii].phi = padLocs[ii].phi + DPHI_AT_SMIN;
      }
      else 
      {
        if (isnan(z)) printf("ERROR: ipadlocs %i = %f",ii,z);
        rawCYL[tt][ii].s = getr_fit_eg6(z, tt);

        rawCYL[tt][ii].phi = getdelphi_fit_eg6(z, tt) + padLocs[ii].phi;
        
//        fprintf(stderr,"pad=%8d z=%4.1f tic=%8d s=%12.4f phi=%12.4f \n",ii,z,tt,rawCYL[tt][ii].s,rawCYL[tt][ii].phi);

        if(floor(rawCYL[tt][ii].s) > GEM_RAD || rawCYL[tt][ii].s < CAT_RAD) num_bad++;

        if(rawCYL[tt][ii].s < 0) 
        {
//            fprintf(stderr,"minTDC not reached pad=%d z=%.1f tic=%d  => %f\n",ii,z,tt,rawCYL[tt][ii].s);
            num_neg++;
        }

        if(rawCYL[tt][ii].phi < 0) rawCYL[tt][ii].phi += 2*PI;
        if(rawCYL[tt][ii].phi > 2*PI) rawCYL[tt][ii].phi -= 2*PI;

        if(rawCYL[tt][ii].phi < 0 || rawCYL[tt][ii].phi > 2*PI) 
        {
          printf("ERROR with drift path calculation in init_paths_using_fit");
    	  exit (1);
        }
      }
    }//end loop over time bins
  }//end loop over pads

  return;
}

float getdelphi_fit_eg6(float z, int TDC) {
  float p0,p1,p2,p3,p4;
  float phi;
/*
  p0 =   0.188644  +0.000637086*z    -1.68395e-05*z*z;
  p1 =  -0.0059482 -8.32225e-05*z    +1.60137e-06*z*z;
  p2 =   0.000444128  +3.754e-06*z   -7.23783e-08*z*z;
  p3 =  -7.83113e-06  -7.41789e-08*z +1.31674e-09*z*z;
  p4 =   5.44596e-08  +5.08519e-10*z -7.66962e-12*z*z;
*/
/*
  // new drift paths from 1.27  and from 1.2 GeV data (from pass1v5):: 02/20/2014
   p0 = 0.174993     +0.000909321*z   +8.13101e-06*z*z   -2.82105e-07*z*z*z;
   p1 =-0.00465778   -0.000104652*z   -9.15424e-07*z*z   +2.83384e-08*z*z*z;
   p2 = 0.00041491   +4.02323e-06*z   +1.87709e-08*z*z   -9.83958e-10*z*z*z;
   p3 =-7.831e-06    -7.18259e-08*z   -1.28834e-10*z*z   +1.54054e-11*z*z*z;
   p4 = 5.7514e-08   +4.6829e-10 *z   +5.66278e-13*z*z   -8.84286e-14*z*z*z;
*/

   // 04/04/2014: new drift paths extracted from pass1v6 after correcting the density unit in the simulation
   p0 =  0.14222       -6.52562e-05*z   +4.06768e-06*z*z; 
   p1 = -0.00147368   +5.64924e-06*z   -7.31944e-07*z*z;
   p2 =  0.000216222   +6.25749e-09*z   +1.8923e-08*z*z;
   p3 = -3.82450e-06  -6.29825e-09*z   -1.89627e-10*z*z;
   p4 =  3.22973e-08   +7.52017e-11*z   +1.08564e-12*z*z;
 
  phi =  p0 + p1*TDC + p2*TDC*TDC + p3*TDC*TDC*TDC+ p4*TDC*TDC*TDC*TDC;


  return phi;
}


/*
float getr_fit_eg6(float z, int TDC)
{
  float DriLen,TotTime,DriSpe,CurR,DPhi,TempR[NAL_SAMP];
  int CurT,jj;
  const float TDCmax = 72.6 - 0.0044*z - 7.8e-04*z*z;
  float maxrad=GEM_RAD;
  if (fabs(z-42.5)<1.e-3 && TDC==15) maxrad += 0.1;
  if (TDC > TDCmax) return 0.; // hit will be rejected by clean_hits
  if (TDC < TPC_TZERO) return 9999.; // hit will be rejected by clean_hits
  for (jj = TPC_TZERO; jj<=TDCmax ; jj++) TempR[jj] = -1;
  
  DriLen = 36; // where did this number come from?
  TotTime = TDCmax - TPC_TZERO;
  if (TotTime <= 0) return -1;
  
    DriSpe = 0.618483 +4.0334e-05*z +8.35178e-06 *z*z;
    CurR = CAT_RAD;
    CurT = TDCmax;
    DriLen = 0.;
    TempR[CurT] = CurR;

    while (CurR < maxrad)
     {
      DPhi = getdelphi_fit_eg6(z,CurT-1) - getdelphi_fit_eg6(z,CurT);
      if(DriSpe*DriSpe-pow(CurR*DPhi,2) > 0)
       {
        CurR += sqrt(DriSpe*DriSpe-CurR*CurR*DPhi*DPhi);
        }
      else
       {
        printf("DriSpe too low!!! at z= %f\n",z);
        break;
        }
  
     CurT -= 1;
     DriLen = DriLen + DriSpe;
     TempR[CurT] = CurR;
   }
  
  return TempR[TDC];
 }
*/  


float getr_fit_eg6(float z, int TDC)
{
  float DriLen,TotTime,DriSpe,CurR,DPhi,TempR[NAL_SAMP];
  int CurT,jj;
  float TDCmax = get_TDCmax_eg6(z);
  float maxrad=GEM_RAD;
  if (fabs(z-42.5)<1.e-3 && TDC==15) maxrad += 0.1;
  if (TDC > TDCmax) return 0.; // hit will be rejected by clean_hits
  if (TDC < TPC_TZERO) return 9999.; // hit will be rejected by clean_hits
  for (jj = TPC_TZERO; jj<=TDCmax ; jj++) TempR[jj] = -1;
  
  DriLen = 36; // where did this number come from? why we are using in some files this value to be 38
  TotTime = TDCmax - TPC_TZERO;
  if (TotTime <= 0) return -1;
  
  for (jj = 0 ; jj<5 ; jj++){
    DriSpe = DriLen / TotTime;
//    DriSpe = 0.60911 +6.04706e-05*z + 6.43316e-06*z*z;
    CurR = CAT_RAD;
    CurT = TDCmax;
    DriLen = 0.;
    if (jj==4) TempR[CurT] = CurR;  
    TempR[CurT] = CurR;
    while (CurR < maxrad){
     DPhi = getdelphi_fit_eg6(z,CurT-1) - getdelphi_fit_eg6(z,CurT);
     if(DriSpe*DriSpe-pow(CurR*DPhi,2) > 0) {
       CurR += sqrt(DriSpe*DriSpe-CurR*CurR*DPhi*DPhi);
        }
      else {
       //printf("DriSpe too low!!! at z= %f\n",z); ! Commented
       break;
        }
                             
       CurT -= 1;
       DriLen = DriLen + DriSpe;

       if (jj==4) TempR[CurT] = CurR;
//        TempR[CurT] = CurR; // in the case of fixed drift speed
      }
     } 
    return TempR[TDC];
   }



// function to find TDCmax depending on run number and z 

float get_TDCmax_eg6(float z)
{
  float TDCmax, p0, p1, p2, p3;
  int runNumber = runNum;
  // the new TDCmax = p0 + p1*exp(p2*(z-p3)^2)
         if (61448 <= runNumber && runNumber <= 61481){
             p0 = 1.14312e+03 -1.75217e-02 *runNumber ;
             p1 = 3.27339e+00 +5.32577e-05 *runNumber;
             p2 =-7.55131e-02 +1.22429e-06 *runNumber;
             p3 = 3.76627e+03 -6.14148e-02 *runNumber;
             TDCmax = p0 + p1* exp(p2*(z-p3)*(z-p3));
            }   
    else if (61483 <= runNumber && runNumber <= 61611){
             p0 = 1.21405e+02 -9.40705e-04 *runNumber ;
             p1 = 3.89644e+00 +6.33813e-05 *runNumber;  
             p2 = 4.09308e-03 -6.97839e-08 *runNumber;
             p3 =-9.04583e+02 +1.45062e-02 *runNumber;
             TDCmax = p0 + p1* exp(p2*(z-p3)*(z-p3));
           }
    else if (61612 <= runNumber && runNumber <= 61646){
             p0 = 1.39733e+03 -2.16496e-02 *runNumber ;
             p1 =-3.07845e+02 +5.10814e-03 *runNumber;
             p2 =-8.23774e-02 +1.33384e-06 *runNumber;
             p3 =-9.05752e+02 +1.44872e-02 *runNumber;
             TDCmax = p0 + p1* exp(p2*(z-p3)*(z-p3));
            } 
    else if (61655 <= runNumber && runNumber <= 61779){ 
             p0 = 1.45093e+02 -1.33438e-03 *runNumber ;
             p1 = 1.63746e+02 -2.54273e-03 *runNumber; 
             p2 =-5.10501e-04 +5.64359e-09 *runNumber;
             p3 = 4.26282e+02 -7.12408e-03 *runNumber;
             TDCmax = p0 + p1* exp(p2*(z-p3)*(z-p3));
           }
    else if (61791 <= runNumber && runNumber <= 61930){ 
             p0 = 2.18243e+02 -2.51495e-03 *runNumber ;
             p1 = 4.92691e+01 -6.90443e-04 *runNumber;
             p2 =-1.11909e-02 +1.78407e-07 *runNumber;
             p3 = 4.26297e+02 -7.12383e-03 *runNumber;
             TDCmax = p0 + p1* exp(p2*(z-p3)*(z-p3));  
            }
    else if (61931 <= runNumber && runNumber <= 61961){ 
             p0 = 2.18152e+02 -2.51641e-03 *runNumber ;
             p1 = 4.92921e+01 -6.90070e-04 *runNumber;
             p2 =-1.11766e-02 +1.78639e-07 *runNumber;
             p3 = 4.23668e+02 -7.16628e-03 *runNumber;
             TDCmax = p0 + p1* exp(p2*(z-p3)*(z-p3));
            }

    //for the runs which are not included in the fit: take the original TDCmax
    else { TDCmax = 53.1422 -0.0317579*z +0.00199476*z*z; }


 return TDCmax;
}










