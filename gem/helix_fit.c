#include "rtpc.h"
#include <stdio.h>
#include <math.h>
/* stuff for calling fortran routines and linking cernlib */
//#define  f2cFortran
/* return values ch2ph and ch2z must be pointer */
// void RWFTHL(int npt, float* rf, float* pf, float* wfi, float* zf, float* wzf,
//	     int iopt, float* vv0, float* ee0, float ch2ph, float ch2z, float* del, float* delz) {
//   rwfthl_ (&npt,rf,pf,wfi,zf,wzf,&iopt,vv0,ee0,&ch2ph,&ch2z,del,delz);
//}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void helix_fit(dummy, tracknum, chainnum)
     int dummy, tracknum, chainnum;

     /* fit a helix to the points on chain #chainnum and store results as
	track #tracknum. If (fit_track_to_beamline), then DO include 
	the beamline (x=0,y=0,z=unconstrained) in the fit */
{

  int jj;
  float xx,yy;
  float my_phi;
  float rf[MAX_HITS_ON_CHAIN];
  float pf[MAX_HITS_ON_CHAIN];
  float wfi[MAX_HITS_ON_CHAIN];
  float zf[MAX_HITS_ON_CHAIN];
  float wzf[MAX_HITS_ON_CHAIN];
  int npt, npch2;
  float vv0[5];
  float ee0[15];
  float ch2ph;
  float ch2z;
  float del[MAX_HITS_ON_CHAIN];
  float delz[MAX_HITS_ON_CHAIN];
  float deviationSum;
  //for testing the c version of the helix fit routine
  int   helix_debug = FALSE;
//  int   use_ccode   = TRUE;
  int h_error;

  cylindrical closest_point;
  float sigma, tot_wt;

  if(chainnum>MAX_NUM_CHAINS) errors[5]++;// this is impossible

  if(helix_debug) printf(" helix_fit: fitting %d hits on chain %d to make track %d\n",
         		num_hits_this_chain[chainnum],chainnum,tracknum);	 
  npt=0;
  tot_wt= 0.0;
  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {
    if(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)
    {
      if(GEMDEBUG) printf("fitted hit status= %d\n",hh_hitlist[chain_hits[chainnum][jj]].status);
      rf[npt] = hh_hitlist[chain_hits[chainnum][jj]].r;
      pf[npt] = hh_hitlist[chain_hits[chainnum][jj]].phi;
      zf[npt] = hh_hitlist[chain_hits[chainnum][jj]].z;
      wfi[npt]=  hh_hitlist[chain_hits[chainnum][jj]].q;
      tot_wt+= wfi[npt];
      wzf[npt]=  hh_hitlist[chain_hits[chainnum][jj]].q;
      npt++;
    }
  }
  if(fit_track_to_beamline)
  {
    /* Fit it first using beam (x,y)= (0,0) to get an approximate z. Then fit 
       again using the (x,y) of the known beamline at that z. */
    rf[npt]= 0.0;
    pf[npt]= 0.0;
    zf[npt]= 0.0;
    wfi[npt]= tot_wt;
    wzf[npt]= 0.0; /* zero weight for Z on the beamline point*/
    npt++;
//    if(!use_ccode) RWFTHL(npt,rf,pf,wfi,zf,wzf,rwfopt,vv0,ee0,ch2ph,ch2z,del,delz);
//    else rwfthc(npt,rf,pf,wfi,zf,wzf,rwfopt,&h_error,vv0,ee0,&ch2ph,&ch2z,del,delz); 
    rwfthc(npt,rf,pf,wfi,zf,wzf,rwfopt,&h_error,vv0,ee0,&ch2ph,&ch2z,del,delz); 

    /* Now set up for second fit, this time using corrected beam (x,y) */
    npt--;
    z2beamxy(vv0[4],&xx,&yy); /*gets (x,y) of beam at z -- in BONUS Coordinates */
    rf[npt]= sqrt(xx*xx + yy*yy);
    pf[npt]= atan2(yy,xx);
    zf[npt]= 0.0;
    wfi[npt]= tot_wt;
    wzf[npt]= 0.0; /* zero weight for Z on the beamline point*/
    npt++;
   }

  if(helix_debug)  
  { 
    printf("helix fitting npt = %d\n", npt); 
    for (jj=0; jj<npt; jj++) 
      printf("(%d,r,phi,z)= (%5.4f, %5.4f, %5.4f,%5.0f)\n",jj,rf[jj],pf[jj],zf[jj],wzf[jj]); 
  } 
  //if(!use_ccode) RWFTHL(npt,rf,pf,wfi,zf,wzf,rwfopt,vv0,ee0,ch2ph,ch2z,del,delz);
  //else           rwfthc(npt,rf,pf,wfi,zf,wzf,rwfopt,&h_error,vv0,ee0,&ch2ph,&ch2z,del,delz); 
  rwfthc(npt,rf,pf,wfi,zf,wzf,rwfopt,&h_error,vv0,ee0,&ch2ph,&ch2z,del,delz); 
 
  dzds[tracknum]= vv0[1];
  helix_error[tracknum] = 0;//h_error;
  r_0[tracknum] = -1.0/vv0[0]; /* radius of curvature. +ve for +ve tracks in BoNuS*/
  phi0[tracknum]= vv0[2]; /* in xy plane, direction of track relative to x axis */
  my_phi= phi0[tracknum]+PI;
  if (vv0[0]<0.0) my_phi+=PI;
  if(my_phi>2.0*PI) my_phi-=2.0*PI;
  dca[tracknum]= vv0[3]; /* dca = distance of closest approach to z-axis */
  x_close[tracknum]= -sin(my_phi)*(-vv0[3]);
  y_close[tracknum]=  cos(my_phi)*(-vv0[3]);
  z_close[tracknum] =  vv0[4];     /* z at point of closest approach */
  x_0[tracknum] = -sin(my_phi)*((-vv0[3])+fabs(1.0/vv0[0]));
  y_0[tracknum] =  cos(my_phi)*((-vv0[3])+fabs(1.0/vv0[0]));
  
  deviationSum = 0.0;
  npch2= 0;
  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {
    if(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)
    {
      hh_hitlist[chain_hits[chainnum][jj]].del = 
                  dca_to_helix(r_0[tracknum],  x_0[tracknum], y_0[tracknum], 
                  dzds[tracknum],
                  x_close[tracknum], y_close[tracknum], 
                  z_close[tracknum],
                  &hh_hitlist[chain_hits[chainnum][jj]], 
                  &closest_point, &sigma);
      //deviationSum += hh_hitlist[chain_hits[chainnum][jj]].del;
      deviationSum += pow(hh_hitlist[chain_hits[chainnum][jj]].res_r/SIGMA_R,2);
      deviationSum += pow(hh_hitlist[chain_hits[chainnum][jj]].res_phi/SIGMA_PHI,2);
      deviationSum += pow(hh_hitlist[chain_hits[chainnum][jj]].res_z/SIGMA_Z,2);
      npch2+=3;
    } 
  } 

  if(fit_track_to_beamline) 
  {
    deviationSum += pow(dca[ntracks]/SIGMA_R , 2);
    npch2++;
  }
  
  if(npch2>4) {chi2[tracknum] = deviationSum/(npch2-4);} 
  else {chi2[tracknum]=9999.9;}
  if(helix_debug) 
  {
    printf("deviationsum = %f\n",deviationSum);
    printf("Helix_fit returns 1/r,dzds,phi0,dca,z_close,x_center,y_center\n");
    printf("%6f %6f %6f %6f %6f %6f %6f\n",
   	vv0[0],vv0[1],vv0[2],vv0[3],vv0[4],x_0[tracknum],y_0[tracknum]); 
  }
  return;
}
