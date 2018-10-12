//#include <stdio.h>
#include <math.h>
#include "rtpc.h"

#define f2cFortran
void NSORTR(float* xbd, int id, int nint, int fwd_bkd) {
  nsortr_ (xbd, &id, &nint, &fwd_bkd);
}

int findwall(int tracknum, int fwd_bkd, float *sig1, float *sig2)
/* By projecting the helix described by 'tracknum' either 
   forward (fwd_bkd=+1) or backward (fwd_bkd=-1), find the
   next two points where the track crosses the RTPC sensitive
   volume boundary. The helix "begins" at the point of closest
   approach to the z-axis, at which point 'sigma' is defined
   to be zero. This is called "the dca point". The function
   value returned is the number of boundaries found in the
   direction fwd_bkd, to a maximum of two. The positions of
   those crossings are returned as sig1 and sig2.
   
   sig1 and  sig2  ARE ALWAYS POSITIVE IN THE FWD_BKD  DIRECTION
   SO ACCOUNT FOR ANY DESIRED SIGN CONVENTION IN  THE CALLING CODE.
   
   Should the dca point be within the RTPC sensitive volume,
   it may be useful to re-call this function with fwd_bkd=-1
   in order to find the point where the track entered the
   chamber. */

/* REVISION HISTORY
   2-FEB-07  hcf repaired a bug in finding intersections of helix with cage boards
   2012 nab repaired bug with untiliazed values assigned to sig1,sig2
*/
{
  float sign_q, radius, rcent, phi_dca, ztest;
  float rlim, cg, sg, gamma, mydzds;
  float sig_u, sig_d, sig_i, sig_o, sig;
  float x1, y1, x2, y2;
  float psi1, psi2, dphi;
  int ii;
#define NPMAX 10
  int nint; /* number of points where this helix intersects any chamber boundary */
  float xbd[NPMAX][4]; /* penetration points */

  //////nate added on 12-05-05 to keep z position of exit point/////
  float r_exit, r_point;
  r_exit = 0.0;
  z_exit = -10000;

  nint= 0;
  
  radius= fabs(r_0[tracknum]);
  mydzds= dzds[tracknum];
  sign_q= radius/r_0[tracknum];
  gamma= atan(mydzds);
  cg= cos(gamma);
  sg= sin(gamma);

  /* rcent = R  at helix axis */
  rcent= sqrt( x_0[tracknum]*x_0[tracknum]  + y_0[tracknum]*y_0[tracknum] );
  /* angle of line  from helix center to dca */
  phi_dca= atan2(y_close[tracknum]-y_0[tracknum],x_close[tracknum]-x_0[tracknum]);

  /* Find boundary where helix enters sensitive volume of rtpc...*/
  /* where it intersects the upstream wall at z=-100: */
  ztest=-ZEND;
  sig_u= fwd_bkd*(ztest-z_close[tracknum])/(mydzds*cg);
  sig= sig_u;
  xbd[nint][1]= radius*cos(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+x_0[tracknum];
  xbd[nint][2]= radius*sin(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+y_0[tracknum];
  xbd[nint][3]= z_close[tracknum] + fwd_bkd*sig*sg;
  xbd[nint][0]= sig;
  if( (sig>0.0) && inside_rtpc(xbd[nint][1],xbd[nint][2],xbd[nint][3])) nint++;

  /* where it intersects the downstream wall at z=+100: */
  ztest=+ZEND;
  sig_d= fwd_bkd*(ztest-z_close[tracknum])/(mydzds*cg);
  sig= sig_d;
  xbd[nint][1]= radius*cos(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+x_0[tracknum];
  xbd[nint][2]= radius*sin(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+y_0[tracknum];
  xbd[nint][3]= z_close[tracknum] + fwd_bkd*sig*sg;
  xbd[nint][0]= sig;
  if( (sig>0.0) && inside_rtpc(xbd[nint][1],xbd[nint][2],xbd[nint][3])) nint++;

  /* whether it intersects the inner radius (cathode) at r=30.0 */
  sig_i= -999.0;
  rlim= CAT_RAD;
  if(circle_int(radius,x_0[tracknum],y_0[tracknum],rlim,0.0,0.0,&x1,&y1,&x2,&y2))
  {
    /* helix's internal angle starts at dca and runs to these two points */
    /* Angles in the lab system are... */
    psi1= atan2(y1-y_0[tracknum], x1-x_0[tracknum]);
    psi2= atan2(y2-y_0[tracknum], x2-x_0[tracknum]);
    /* Move angles to the helix system */
    for(psi1-= phi_dca;psi1<0.0;psi1+=2.0*PI);
    for(psi2-= phi_dca;psi2<0.0;psi2+=2.0*PI);
    /*  If the particle is +ve, choose  +ve turning angles, etc. */
    if(fwd_bkd*sign_q < 0.0) 
    {
      psi1= 2.0*PI  - psi1;
      psi2= 2.0*PI  - psi2;
    }

    /* We care only about the first penetration */
    dphi= psi1; 
    if (psi2 < dphi) dphi= psi2;

    sig= fabs(radius*dphi/cg);
    xbd[nint][1]= radius*cos(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+x_0[tracknum];
    xbd[nint][2]= radius*sin(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+y_0[tracknum];
    xbd[nint][3]= z_close[tracknum] + fwd_bkd*sig*sg;
    xbd[nint][0]= sig;
    if(inside_rtpc(xbd[nint][1],xbd[nint][2],xbd[nint][3])) nint++;
    sig_i= sig;
  }
 
  /* whether it intersects the outer radius (GEM1) at r=60.0 */
  sig_o= -999.0;
  rlim= GEM_RAD;
  if(circle_int(radius,x_0[tracknum],y_0[tracknum],rlim,0.0,0.0,&x1,&y1,&x2,&y2))
  {
    /* helix's internal angle starts at dca and runs to these two points */
    /* Angles in the lab system are... */
    psi1= atan2(y1-y_0[tracknum], x1-x_0[tracknum]);
    psi2= atan2(y2-y_0[tracknum], x2-x_0[tracknum]);
    /* Move angles to the helix system */
    for(psi1-= phi_dca;psi1<0.0;psi1+=2.0*PI);
    for(psi2-= phi_dca;psi2<0.0;psi2+=2.0*PI);
    /*  If the particle is +ve, choose  +ve turning angles, etc. */
    if(fwd_bkd*sign_q < 0.0) 
    {
      /* a little thought will prove that this is redundant, I think. BUt for now... */
      psi1= 2.0*PI - psi1;
      psi2= 2.0*PI - psi2;
    }
    /* We care only about the first penetration */
    dphi= psi1; if (psi2 < dphi) dphi= psi2;

    sig= fabs(radius*dphi/cg);
    xbd[nint][1]= radius*cos(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+x_0[tracknum];
    xbd[nint][2]= radius*sin(fwd_bkd*sign_q*sig*cg/radius+phi_dca)+y_0[tracknum];
    xbd[nint][3]= z_close[tracknum] + fwd_bkd*sig*sg;
    xbd[nint][0]= sig;
    if(inside_rtpc(xbd[nint][1],xbd[nint][2],xbd[nint][3])) nint++;

    sig_o= sig;
  }

  /* order the points by path length */
  if(nint>1) NSORTR(xbd,4,nint,fwd_bkd);

  for (ii=0; ii<nint; ii++)
  {
    r_point = sqrt(xbd[ii][1]*xbd[ii][1] + xbd[ii][2]*xbd[ii][2]);
    if(r_point > r_exit)
    {
      r_exit = r_point;
      z_exit = xbd[ii][3];
    }
  }
  
  *sig1 = nint>0 ? xbd[0][0] : -9999;
  *sig2 = nint>1 ? xbd[1][0] : -9999;
 
  return (nint);
}
