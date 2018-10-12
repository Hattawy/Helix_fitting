#include "rtpc.h"
#include <stdio.h>

// nab modified 2012

int countPads(int tracknum, int chainnum)
{
  /////////////////////////////////////////////////////////////////////////
  //this is a simple routine to count how on many channels the track induced a
  //signal
  ////////////////////////////////////////////////////////////////////////

  int padList[NUM_PADS];
  int numPads, foundPad, thisPad;
  int jj, kk;

  numPads = 0;
  foundPad = FALSE;

  if(chainnum>MAX_NUM_CHAINS) printf("chainnum is to big in countPads()!\n");

  //initialize the pad list
  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {

    if(!(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)) continue;

    thisPad = hh_hitlist[chain_hits[chainnum][jj]].pad;

    foundPad = FALSE;
   
    for (kk=0; kk<numPads; kk++)
    {
      if (padList[kk]==thisPad) 
      {
        foundPad=TRUE;
        break;
      }
    }
    
    if(!foundPad) 
    {
      padList[numPads] = thisPad;
      numPads++;
    }

  }//end of loop over hits in the chain

  return (numPads);
}
