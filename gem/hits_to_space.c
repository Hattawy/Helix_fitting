//#include <stdio.h>
#include "rtpc.h"

/* Hits on hitlist have only fundamental data (pad number and altro time)
   stored. Add in the spatial information. */
void hits_to_space()
{
  int ii, tic, pad_num;
  for( ii=0; ii<hh_num_hits; ii++)
  {
   /* printf("%d %d %d \n",ii ,hh_num_hits, HH_MAX_NUM_HITS);
    printf("%5f %5f %5f %5f %5f \n", tic, pad_num, hh_hitlist[ii].z, hh_hitlist[ii].phi, hh_hitlist[ii].r);*/
    
    //if (hh_num_hits>= HH_MAX_NUM_HITS) break; (nab commented out 2012)

    tic=  hh_hitlist[ii].t;
    pad_num= hh_hitlist[ii].pad;
    hh_hitlist[ii].z= rawCYL[tic][pad_num].z;
    hh_hitlist[ii].r= rawCYL[tic][pad_num].s;
    hh_hitlist[ii].phi= rawCYL[tic][pad_num].phi;
   /* if ( rawCYL[tic][pad_num].s <= 0 )
    {
      printf("%d %d %d \n",ii ,hh_num_hits, HH_MAX_NUM_HITS);
      printf("%5f %5f %5f %5f %5f \n", tic, pad_num, hh_hitlist[ii].z, hh_hitlist[ii].phi, hh_hitlist[ii].r);
   }*/
  }
  return;
}
   
