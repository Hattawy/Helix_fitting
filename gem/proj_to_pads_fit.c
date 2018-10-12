#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rtpc.h"
//
// nab modified 2012
// - tic initialization and convergence loop fixed
//
// return value:
// 0 == good result
// 1 == bad result
//
int proj_to_pads_fit_eg6(cylindrical point, cylindrical* atpad, float *altro_time_at_pad)
{
  float del_phi, myphi, delta, mindelta;
  int tic, niter, flipped, lastmove, mintic;
  const float gooddelta=0.4;

  if (isnan(point.z)) return 1;


// This convergence loop had a few problems:
/*  
  tic = 15+ (int) (point.s-30)/30*55;
  delta =.5;
  niter = 0;
  while(abs(delta)>gooddelta)
  {
    niter++;
    delta = getr_fit_eg6(point.z, tic) - point.s;
    if      (delta<=.4) tic -= 1;
    else if (delta>=.4) tic += 1;
    if (tic<=13 || tic >=75 || niter>30 ) break;
  }
*/

// Replacing with this one:

  niter=0;
  flipped=0;
  lastmove=0;
  mindelta=9999.;
  mintic=-1;

  // initial guess:
  tic = TPC_TZERO + (int)
    ( (GEM_RAD-point.s) / (GEM_RAD-CAT_RAD) * (TPC_TCAT-TPC_TZERO) );

  delta = getr_fit_eg6(point.z, tic) - point.s;

  while ( fabs(delta) > 0.001 )
  {

    if (delta < 0) {
      if (lastmove > 0) flipped++;
      tic -= 1;
      lastmove = -1;
    }
    else
    {
      if (lastmove < 0) flipped++;
      tic += 1;
      lastmove = 1;
    }

    delta = getr_fit_eg6(point.z, tic) - point.s;

    if (fabs(delta) < mindelta) 
    {
      mintic   = tic;
      mindelta = fabs(delta);
    }

    if (flipped > 1) // don't reverse direction more than once
    {
      if (mindelta <= gooddelta) break;
      return 1;
    }
    
    if (niter > 2 && (tic <= TPC_TZERO-2 || tic >= TPC_TMAX)) return 1;

    if (niter > 60) return 1;

    niter++;

  }


  *altro_time_at_pad = mintic;
  
  del_phi = getdelphi_fit_eg6(point.z,mintic);

  
  if(isnan(del_phi) != 0 || del_phi < -99)
  {
    if(GEMDEBUG) printf("del_phi is crazy in proj_to_pads()\n"); 
    return 1;
  }

  myphi= point.phi - del_phi;

  if      (myphi <  0.0)    myphi += 2.0*PI;
  else if (myphi >= 2.0*PI) myphi -= 2.0*PI;
  
  (*atpad).s   = PAD_RAD;
  (*atpad).phi = myphi;
  (*atpad).z   = point.z;

  return 0;
}
