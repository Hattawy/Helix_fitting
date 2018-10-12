//#include <stdio.h>
#include "rtpc.h"
//
// nab modified 2012
// - changed radius tests to be <= or >= instead of ==
//   -- this goes with init_paths_using_fit_eg6 
//
int clean_hits()
     /* Here is where we disable hits for various reasons
	the first of which is that a hit is out of time. */
{
     int ii,jj;
     int nhpp[3200];

     for (ii=0; ii<3200 ; ii++)
       {
       nhpp[ii] = 0;
       }
     
     for (ii=0; ii<hh_num_hits; ii++)
     {
       /* Reset status retaining only HSMALLQ bit (which might have been set in data-input routine */
       hh_hitlist[ii].status &= HSMALLQ;
       
       /* the following two seem redundant but are not because the possible ePath length varies with z */
       if (hh_hitlist[ii].t < TPC_TMIN) hh_hitlist[ii].status |= HDISEAR;
       if (hh_hitlist[ii].t > TPC_TMAX) hh_hitlist[ii].status |= HDISLAT;
       if (hh_hitlist[ii].q < THRESHA) hh_hitlist[ii].status |= HDISEAR;

       // There are cases where drift paths calibration returns radius
       // slightly larger than 1st gem radius.  We are allowing this here:  (nab 2012)

       //if (hh_hitlist[ii].r == SMAX) hh_hitlist[ii].status |= HDISEAR; /* can be set to SMAX in read_path.c */
       //if (hh_hitlist[ii].r >= SMAX) hh_hitlist[ii].status |= HDISEAR;
       if (hh_hitlist[ii].r >= GEM_RAD+0.1) hh_hitlist[ii].status |= HDISEAR;
       
       //if (hh_hitlist[ii].r == SMIN) hh_hitlist[ii].status |= HDISLAT; /* can be set to SMIN in read_path.c */       
       //if (hh_hitlist[ii].r <= SMIN) hh_hitlist[ii].status |= HDISLAT;
       if (hh_hitlist[ii].r <= CAT_RAD-0.1) hh_hitlist[ii].status |= HDISLAT;
       
       /* disable any hits occuring on pads in the hotPads list */
       if (dead_pad(hh_hitlist[ii].pad)) hh_hitlist[ii].status |= HBADPAD;
       
       nhpp[hh_hitlist[ii].pad] += 1;
     }

     /* disable any pad that fire more than HITTHR times in a single event*/
     for (ii=0; ii<3200 ; ii++)
     {
       if (nhpp[ii] <= HITTHR) continue;
       for(jj=0; jj<hh_num_hits; jj++)
       {
         if(hh_hitlist[jj].pad==ii) hh_hitlist[jj].status |= HBADPAD;
       }
     }

     // This will only reject hits by setting their status
     RejectNoise();

     // This can either reject hits, or perform a pedestal correction.
     // MUST be called after scale_hits, because it will reset hh_hitlist.q
     // in case of pedestal subtraction:
     BusyBoards();

     return 1;
}
