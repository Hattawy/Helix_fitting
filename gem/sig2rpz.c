//#include <stdio.h>
#include <math.h>
#include "rtpc.h"

/* function to calculate (r,phi,z) of the point on the helix already
fitted to track 'tracknum' at position 'sigma' along the helix. As
in other routines, sigma is defined to start at the point of closest
approach to the z-axis, given by z_close[tracknum]. The resulting
space point (r,phi,z) is returned in 'point' */

void sig2rpz(int tracknum, float sigma, cylindrical* point)
{ 
  float radius, sign_q, gamma, cg, sg, phi_dca, xx, yy;
  radius= fabs(r_0[tracknum]);
  sign_q= radius/r_0[tracknum];
  gamma= atan(dzds[tracknum]);
  cg= cos(gamma);
  sg= sin(gamma);
  /* angle of line  from helix center to dca */
  phi_dca= atan2(y_close[tracknum]-y_0[tracknum],x_close[tracknum]-x_0[tracknum]);
      
  xx= radius*cos(sign_q*sigma*cg/radius+phi_dca)+x_0[tracknum];
  yy= radius*sin(sign_q*sigma*cg/radius+phi_dca)+y_0[tracknum];

  (*point).s= sqrt(xx*xx + yy*yy);
  (*point).phi= atan2(yy,xx);
  (*point).z= z_close[tracknum] + sigma*sg;
  return;  
}
