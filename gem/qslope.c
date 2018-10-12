#include "rtpc.h"
#include <stdio.h>
#include <math.h>

//
// Calculate ADC-slope for one track (new in 2014, NAB)
//
// INPUT:
// smin/smax     = distance from vertex along helix of first/last hits
// npt           = number of hits on the track
// indices[npts] = index of hits in hh_hitlist array 
//
// OUTPUT:  (via global trackv[] members)
// alpha = slope 
// Rxy   = linear correlation parameter 
// chisq = chi-square (incorrect, but maybe still useful)
// qint  = total track charge from integrating line
//
// NOTE:  This may modify the input array "indices" by sorting it
//

//#define __NPATHBINS__ 12
#define __NPATHBINS__ 7
#define __NMINBINS__ 3

//#define __SLOPFILE__ 1
#define __SLOPFILE__ 0

void qslope(const int tracknum,const float smin,const float smax,
            const int npt,int indices[MAX_HITS_ON_CHAIN],
            const float adcrejfrac)
{
#if __SLOPFILE__
  static FILE* slopfile=NULL;
  if (!slopfile) slopfile=fopen("slopfile.txt","w");
  fprintf(slopfile,"AAAAA %d\n",clas_global_eid);
#endif
  
  int ii;
  float dist;
  int bins[MAX_HITS_ON_CHAIN];
  int nhitsbin[__NPATHBINS__]={};
  float qbin[__NPATHBINS__]={};
  float distbin[__NPATHBINS__]={};

  // initialize globals:
  trackv[tracknum].Rxy = -999999;
  trackv[tracknum].alpha = -999999;
  trackv[tracknum].beta = -999999;
  trackv[tracknum].chisq = -999999;
  trackv[tracknum].qint = 0;
  for (ii=0; ii<npt; ii++)
  {
      hh_hitlist[indices[ii]].res_q = -999999;
      hh_hitlist[indices[ii]].res_qr = -999999;
  }

  const float binsize = (smax-smin) / __NPATHBINS__;

  // number of hits to ignore:
  int nignore = adcrejfrac*npt;
  if (adcrejfrac > 1e-5 && nignore == 0) nignore=1;

  // don't do anything without sufficient data:
  if (npt-nignore < 2) return;

  // if we're truncating, sort indices by charge:
  if (adcrejfrac > 1e-5)
  {
      int swapped=1;
      while (swapped)
      {
          swapped=0;
          for (ii=1; ii<npt; ii++)
          {
              if (hh_hitlist[indices[ii]].q <
                  hh_hitlist[indices[ii-1]].q)
              {
                  const int itmp = indices[ii];
                  indices[ii] = indices[ii-1];
                  indices[ii-1] = itmp;
                  swapped = 1;
              }
          }
      }
  }

  // sum charge in bins of path length:
  for (ii=0; ii<npt-nignore; ii++)
  {
      dist = hh_hitlist[indices[ii]].dist;

      if (fabs(dist-smax)<1e-5) bins[ii] = __NPATHBINS__-1;
      else                      bins[ii] = (dist-smin)/binsize;

      if (bins[ii]<0 || bins[ii]>=__NPATHBINS__)
      {
          fprintf(stderr,"SUMTRACKQ ERROR:  badbin %f %f %f\n",smin,smax,dist);
      }
      else
      {
          nhitsbin[bins[ii]]++;
          qbin[bins[ii]] += hh_hitlist[indices[ii]].q;
          distbin[bins[ii]] += dist;

#if __SLOPFILE__
          fprintf(slopfile,"%d %d %f %f\n",bins[ii],nhitsbin[bins[ii]],
                  hh_hitlist[indices[ii]].q,hh_hitlist[indices[ii]].dist);
#endif
      }
  }
  
  // count filled bins: 
  int nfilledbins=0;
  for (ii=0; ii<__NPATHBINS__; ii++)
      if (nhitsbin[ii]>0) nfilledbins++;
  if (nfilledbins < __NMINBINS__) return;

  // get averages for each bin:
  for (ii=0; ii<__NPATHBINS__; ii++) 
  {
      if (nhitsbin[ii]<=0) continue;
      distbin[ii] /= nhitsbin[ii];
//      qbin[ii] /= nhitsbin[ii]; // use average charge in each bin
  }

/*
  // get std deviations:
  float edistbin[__NPATHBINS__]={};
  for (ii=0; ii<npt-nignore; ii++)
  {
      if (bins[ii]<0 || bins[ii]>=__NPATHBINS__) continue;
      dist = hh_hitlist[indices[ii]].dist;
      edistbin[bins[ii]] += pow(dist - distbin[bins[ii]] ,2);
  }
  for (ii=0; ii<__NPATHBINS__; ii++) 
  {
      if (nhitsbin[ii]<=0) continue;
      edistbin[ii] = sqrt( edistbin[ii] / (nhitsbin[ii]-1) );
  }
*/

#if __SLOPFILE__
  for (ii=0; ii<__NPATHBINS__; ii++)
  {
      fprintf(slopfile,"%5d %12.3f %12.3f\n",nhitsbin[ii],distbin[ii],qbin[ii]);
  }
#endif




  // Simple linear regression:

  // compute means:
  float xmean=0,ymean=0,x2mean=0,y2mean=0,xymean=0;
  for (ii=0; ii<__NPATHBINS__; ii++)
  {
      if (nhitsbin[ii]<=0) continue;
      ymean += qbin[ii];
      xmean += distbin[ii];
      y2mean += pow(qbin[ii],2);
      x2mean += pow(distbin[ii],2);
      xymean += qbin[ii]*distbin[ii];
  }
  xmean /= nfilledbins; 
  ymean /= nfilledbins;
  x2mean /= nfilledbins; 
  y2mean /= nfilledbins;
  xymean /= nfilledbins;

  // compute standard deviations:
  float ystddev=0,xstddev=0;
  for (ii=0; ii<__NPATHBINS__; ii++)
  {
      if (nhitsbin[ii]<=0) continue;
      xstddev += pow(distbin[ii]-xmean,2);
      ystddev += pow(qbin[ii]-ymean,2);
  }
  xstddev = sqrt( xstddev / (nfilledbins-1) );
  ystddev = sqrt( ystddev / (nfilledbins-1) );


  // degree of linear correlation, in range (-1,1):
  const float Rxy = (xymean - xmean*ymean) / 
      sqrt( (x2mean-xmean*xmean) * (y2mean-ymean*ymean) );

  // linear fit parameters:
  const float beta = Rxy * ystddev / xstddev; // slope
  const float alpha = ymean - beta * xmean; // y-intercept

  // calculate chi-square
  float linreg,chisq=0;
  for (ii=0; ii<__NPATHBINS__; ii++)
  {
      if (nhitsbin[ii]<=0) continue;
      linreg = alpha + beta * distbin[ii];
      chisq += pow(qbin[ii]-linreg,2) / fabs(linreg);
  }
  chisq /= (nfilledbins-2);



  // save result of linear regression in globals:
  trackv[tracknum].Rxy = Rxy;
  trackv[tracknum].alpha = alpha;
  trackv[tracknum].beta = beta;
  trackv[tracknum].chisq = chisq;
  // integral from smin to smax:
  trackv[tracknum].qint = alpha*(smax-smin) + beta/2*(smax*smax-smin*smin); 
  // divide integral by the bin width:
  trackv[tracknum].qint /= binsize;





#if __SLOPFILE__
  fprintf(slopfile,"%12.4f %12.4f %12.4f %12.4f\n",alpha,beta,Rxy,chisq);
#endif




  // compute charge-residuals for each pad in the track:
  int pad;
  float expected;
  int   hits[NUM_PADS]={};
  float chrg[NUM_PADS]={};
  float dsum[NUM_PADS]={};
  float dmin[NUM_PADS];
  float dmax[NUM_PADS];
  for (ii=0; ii<NUM_PADS; ii++)
  {
      dmin[ii] =  99999;
      dmax[ii] = -99999;
  }
  for (ii=0; ii<npt; ii++)
  {
      pad = hh_hitlist[indices[ii]].pad;
      if (pad<0 || pad>=NUM_PADS) continue;
      chrg[pad] += hh_hitlist[indices[ii]].q;
      dsum[pad] += hh_hitlist[indices[ii]].dist;
      hits[pad]++;
      if (hh_hitlist[indices[ii]].dist < dmin[pad])
          dmin[pad] = hh_hitlist[indices[ii]].dist;
      if (hh_hitlist[indices[ii]].dist > dmax[pad])
          dmax[pad] = hh_hitlist[indices[ii]].dist;
  }
  for (ii=0; ii<npt; ii++)
  {
      pad = hh_hitlist[indices[ii]].pad;
      if (dmax[pad]-dmin[pad]>0 && hits[pad]>0 && binsize>0)
      {
          expected = (alpha + beta * (dsum[pad]/hits[pad]) ) / binsize;

          hh_hitlist[indices[ii]].res_q =
              chrg[pad] / (dmax[pad]-dmin[pad]) - expected;

          hh_hitlist[indices[ii]].res_qr =
              (chrg[pad] / (dmax[pad]-dmin[pad])) / expected;
      }
  }




}


