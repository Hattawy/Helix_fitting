#include "rtpc.h"
#include <stdio.h>
#include <math.h>

//  June 2014: Adding calls for additional energy calculations (NAB)

float sumtrackq(int tracknum, float *smin, float *smax, int *t_smin, int *t_smax)
     /* returns  function  value as sum of charges in  hits  assigned  to tracknum.
      smin/smax  are returned set to  path  length from  helix  starting point to
      nearest and  farthest measured points,  respectively 
     t_smin/t_smax are the corresponding altro-times.*/
{
  int jj,npt;
  cylindrical closest_point;
  float sigma,sumq;
  int indices[MAX_HITS_ON_CHAIN];

  const int chainnum = track_chain[tracknum];

  npt=0;
  sumq= 0.0;
  *smin= 9999999.9;
  *smax=-9999999.9;
 
  if(GEMDEBUG)   printf("start sumtrackq()...\n"); 
  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {
    if(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)
    {
      sumq += hh_hitlist[chain_hits[chainnum][jj]].q;
      dca_to_helix(r_0[tracknum],  x_0[tracknum], y_0[tracknum], dzds[tracknum],
    	       x_close[tracknum], y_close[tracknum], z_close[tracknum],
    	       &hh_hitlist[chain_hits[chainnum][jj]], &closest_point, &sigma);

      if(sigma < *smin)
      {
        *smin= sigma;
        tmax  = hh_hitlist[chain_hits[chainnum][jj]].t;
        *t_smin= tmax;
      }
      if(sigma > *smax)
      {
        *smax= sigma;
        tmin  = hh_hitlist[chain_hits[chainnum][jj]].t;
        *t_smax= tmin;
      }
      
      // store info for upcoming calculations:
      hh_hitlist[chain_hits[chainnum][jj]].dist = sigma;
      indices[npt++] = chain_hits[chainnum][jj];
    }
  }

  
  // Advanced energy calculations:

  // truncated dedxs (NOTE, this reorders "indices"):
  const float adcfrac = qtrunc(tracknum,*smin,*smax,npt,indices);

  // calculate slope of ADC vs path length:
  qslope(tracknum,*smin,*smax,npt,indices,0.);

//  // ADC-slope with ADC-truncation:
//  qslope(tracknum,*smin,*smax,npt,indices,adcfrac);



  if(GEMDEBUG)   printf("finish sumtrackq()...\n"); 
  return(sumq);

}
  
