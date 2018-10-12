#include "rtpc.h"
#include <stdio.h>

//
// Calculate a few truncated dedx numerators (new in 2014, NAB)
//
// INPUT:
// smin/smax     = distance from vertex along helix of first/last hits
// npt           = number of hits on the track
// indices[npts] = index of hits in hh_hitlist array
//
// OUTPUT:  (via global trackv[] members)
// qtruncA  = total track charge ignoring large ADCs
// qtruncL  = total track charge ignoring ends of track length
// qtruncAL = total track charge ignoring both
//
// NOTE:  This may modify the input array by sorting it
//

float qtrunc(const int tracknum,const float smin,const float smax,
             const int npt,int indices[MAX_HITS_ON_CHAIN])
{
  // initialize globals:
  trackv[tracknum].qtruncA=0;
  trackv[tracknum].qtruncL=0;
  trackv[tracknum].qtruncAL=0;
  trackv[tracknum].sminL=-99999;
  trackv[tracknum].smaxL=-99999;

  // fractions of hits to ignore for truncations:
  const float adcfrac=0.30;  // highest-ADC hits
  const float lenfrac=0.10;  // hits at ends of track
  const float adcfrac2=0.20; // highest-ADC hits when combined with lenfrac

  int jj,itmp,swapped,nignore,nignore2;

  if (npt < 4) return 0;

  //////////////// ADC-truncated track charge ////////////////////////////
  {
      // bubble-sort hits by increasing charge:
      swapped=1;
      while (swapped)
      {
          swapped=0;
          for (jj=1; jj<npt; jj++)
          {
              if (hh_hitlist[indices[jj]].q < 
                  hh_hitlist[indices[jj-1]].q)
              {
                  itmp = indices[jj];
                  indices[jj] = indices[jj-1];
                  indices[jj-1] = itmp;
                  swapped = 1;
              }
          }
      }
      // calculated track charge, ignoring hits with highest charge:
      nignore = adcfrac * npt;
      if (nignore==0) nignore=1;
      if (npt-nignore > 1)
      {
          for (jj=0; jj<npt-nignore; jj++)
          {
              trackv[tracknum].qtruncA += hh_hitlist[indices[jj]].q;
              hh_hitlist[indices[jj]].status |= HTRUNCA;
          }
      }
  }


  ///////////// Path-length-truncated track charge ////////////////////////
  {
      // bubble-sort hits by increasing distance along track:
      swapped=1;
      while (swapped)
      {
          swapped=0;
          for (jj=1; jj<npt; jj++)
          {
              if (hh_hitlist[indices[jj]].dist < 
                  hh_hitlist[indices[jj-1]].dist)
              {
                  itmp = indices[jj];
                  indices[jj] = indices[jj-1];
                  indices[jj-1] = itmp;
                  swapped = 1;
              }
          }
      }
      // calculated track charge, ignoring hits at ends of track:
      nignore = lenfrac * npt;
      if (nignore==0) nignore=1;
      if (npt-2*nignore > 1)
      {
          for (jj=nignore; jj<npt-nignore; jj++)
          {
              trackv[tracknum].qtruncL += hh_hitlist[indices[jj]].q;
              hh_hitlist[indices[jj]].status |= HTRUNCL;
          }

          // Save new smin and smax after path-length truncation:
          trackv[tracknum].sminL = hh_hitlist[indices[nignore]].dist;
          trackv[tracknum].smaxL = hh_hitlist[indices[npt-nignore-1]].dist;
      }
/*
      // calculated track charge, ignoring ends of track based on path length fraction:
      const float smin2 = smin + (smax-smin)*lenfrac;
      const float smax2 = smax - (smax-smin)*lenfrac;
      for (jj=0; jj<npt; jj++) 
      {
          if (distance[jj] < smin2) continue;
          if (distance[jj] > smax2) break;
          trackv[tracknum].qtruncL2 += hh_hitlist[indices[jj]].q;
      }
*/
  }

  //////////// ADC- and Path-Length- truncated track charge /////////////////
  {
      // # of hits to ignore at beginning and end of track:
      nignore = lenfrac * npt;
      if (nignore==0) nignore=1;
      
      // # of hits to ignore for ADC-truncataion:
      nignore2 = adcfrac2 * (npt-nignore*2);
      if (nignore2==0) nignore2=1;

      if (npt-2*nignore-nignore2 > 0)
      {
          // NOTE, this relies on already being sorted by distance.
          // Now, sort by increasing ADC only in middle of track:
          swapped=1;
          while (swapped)
          {
              swapped=0;
              for (jj=nignore+1; jj<npt-nignore; jj++)
              {
                  if (hh_hitlist[indices[jj]].q < 
                      hh_hitlist[indices[jj-1]].q)
                  {
                      itmp = indices[jj];
                      indices[jj] = indices[jj-1];
                      indices[jj-1] = itmp;
                      swapped = 1;
                  }
              }
          }
          // calculate track charge, ignore ends of track and high-ADCs:
          for (jj=nignore; jj<npt-nignore-nignore2; jj++)
          {
              trackv[tracknum].qtruncAL += hh_hitlist[indices[jj]].q;
              hh_hitlist[indices[jj]].status |= HTRUNCAL;
          }
      }
  }


  return adcfrac;

}

