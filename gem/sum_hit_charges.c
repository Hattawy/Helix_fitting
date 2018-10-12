#include <stdio.h>
#include "rtpc.h"

// nab modified 2012:
// - removed unused variables & stopped passing globals

void sum_hit_charges(int *npads_hit, float *qtot)
{
  int jj,pad_num,time;
  float charge;
  float prow, pcol;

  *npads_hit= 0;
  *qtot= 0.0;
      
  for (jj=0; jj< hh_num_hits; jj++)
  {
//   if ( ! (hh_hitlist[jj].status & HISUSED) ) continue;

    pad_num= hh_hitlist[jj].pad;
    charge= hh_hitlist[jj].q;
    time= hh_hitlist[jj].t;

    *qtot += charge;

    if      (time < TPC_TGEM) qtot_early += charge;
    else if (time > TPC_TCAT) qtot_late  += charge;	  
    else                      qtot_intime+= charge;

    if((charge > THRESHA) && (hh_hitlist[jj].status & HISUSED)) *npads_hit += 1;

    if(pad_num != prowcol(padLocs[pad_num].z,padLocs[pad_num].phi,&prow,&pcol))
    {
      errors[12]++;
      
      printf(">>> pad_num disagreement in sum_hit_charges! This should not happen. \n");
      printf(">>> pad %d is at (%f,%f) but prowcol gives pad %d\n",
    	 pad_num,padLocs[pad_num].z,padLocs[pad_num].phi,
    	 prowcol(padLocs[pad_num].z,padLocs[pad_num].phi,&prow,&pcol));
    }
  }
  
  return;
}
