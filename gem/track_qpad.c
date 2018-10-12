#include <stdio.h>
#include "rtpc.h"

void track_qpad(int tracknum, float pp, float dqdx, float sdist, float vtl)
     /* for each pad contributing to the track, computes and fills 
	histogram with the ratio of that pad's measured charge to
        the average charge/pad along the whole track 
        and fill histogram, for each pad on track, with ratio of
        this track's *dqdx* normalized to expected dqdx given track
        momentum *pp* 
	
	3/17/06 hcf
	Add code to figure out the track length subtended by each
	pad and make histogram of measured charge / expected charge
	for each pad and track, where expected charge is the total
	charge *(subtended track length / total track visible track
	length). Had to add parameters sdist and vtl, which are the
	starting distance along the helix and the length of the helix
	we wish to examine.

	7/01/07 nkb
	altered to only keep the counter for the number of pad uses on a track
     */
{
  int ii, jj, pad_num, old_pad, time;
  int num_pads, used, num_hits_used;
  int chainnum;
  float sumq, thisq;
  float qtime[NAL_SAMP];

  typedef struct {
    int pad;
    float qsum;
  } listt;
#define MAXPADLIST 1000
  listt pad_list[MAXPADLIST];
  float qlist[MAXPADLIST];

  chainnum= track_chain[tracknum];

  /* Here we're going to step through all hits on the track and sum up
     the charges that land on each pad */
  num_pads=0;
  sumq= 0.0;
  num_hits_used= 0;

  for (ii=0; ii<NAL_SAMP; ii++) {qtime[ii]= 0.0;}
  
  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {
    if(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)
    {
      /* get the total charge on this track */
      thisq= hh_hitlist[chain_hits[chainnum][jj]].q;
      sumq+= thisq;
      if(num_hits_used<MAXPADLIST)
        qlist[num_hits_used]= thisq;
      else
        printf(">>ERROR! Ran out of space for qlist in track_qpad.c\n");

      num_hits_used++;
      
      /* identify this pad */
      pad_num= hh_hitlist[chain_hits[chainnum][jj]].pad;
      time=  hh_hitlist[chain_hits[chainnum][jj]].t;

      qtime[time]+= thisq;

      /* if this pad has already been seen on this track add this hit's
         charge to the sum for this pad */
      used= FALSE;
      for (ii=0; (ii<num_pads)&&(!used); ii++)
      {
        old_pad= pad_list[ii].pad;
        if(old_pad == pad_num)
        {
          pad_list[ii].qsum += thisq;
          used= TRUE;
        }
      }
      /* if it wasn't already seen, create and initialize a new entry */
      if(!used)
      {
        if(num_pads < MAXPADLIST)
    	{
    	  pad_list[num_pads].pad= pad_num;
    	  pad_list[num_pads].qsum = thisq;
    	  num_pads++;
    	}
        else
    	{
    	  printf("Too many pads on padlist in track_qpad. Pad ignored.\n");
    	}
      }
    }
  }
  return;
}
  
