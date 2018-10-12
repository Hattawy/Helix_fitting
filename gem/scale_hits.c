#include <stdio.h>
#include "rtpc.h"

void scale_hits()
{
  int ii,pad;

  for (ii=0; ii<hh_num_hits; ii++)
  {

    pad = hh_hitlist[ii].pad;
    if( pad<0 || pad>NUM_PADS-1)
    {
      errors[11]++;
      continue;
    }
   
    // apply pad-by-pad averaged pedestal
//    hh_hitlist[ii].q -= eg6peds[pad];
//    if (hh_hitlist[ii].q < 1e-6) 
//    {
//      hh_hitlist[ii].status |= HSMALLQ; 
//      continue;
//    }

    // add the amount of signal suppressed by the altro offset
    hh_hitlist[ii].q += altro_offset;
    
    // scale each hit's pulse-height by the gain of the corresponding pad
    hh_hitlist[ii].q /= corr_fac[pad];
  }
  return;
}
