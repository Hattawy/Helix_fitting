//#include <stdio.h>
#include <math.h>
//
// nab modified 2012
// - replaced fortran to get unit vector with just c code
//
float lineandhelix(float xc,float yc, float zh, float rho, float dzds,
		  float x0, float y0, float z0, float dxdz,  float dydz)
     /* given a helical track described by 
	    xc, yc, rho, dzds = center and radius and dzdx of a helix // to z-axis,
	    zh is Z of the point on the helix closest to the Z-axis;
	and a line passing through (x0,y0,z0) with slopes dxdz, dydz,

	return 'sigma'  (distance along helix) of point of closest approach
	of  the helix to the line.  */
{

  float sigma,myrho,sign,mag;
  typedef  float vec[6];  
  vec lin;
  float x0h,y0h,z0h; /* line intercepts in internal helix coordinate system */
  float axh,ayh,azh; /* line direction cosines in internal helix CS */
  float psi, cpsi, spsi, gamma, cg, sg;

  myrho= fabs(rho);
  sign= rho/myrho;

  /* convert line direction parameters into direction  cosines */
  lin[0]= x0;
  lin[1]= y0;
  lin[2]= z0;
  lin[3]= dxdz;
  lin[4]= dydz;
  lin[5]= 1.0;
  /* convert line direction parameters into direction  cosines */
  mag=sqrt( lin[3]*lin[3] + lin[4]*lin[4] + lin[5]*lin[5] );
  lin[3] /= mag;
  lin[4] /= mag;
  lin[5] /= mag;

  /* CS rotation from  universe to helix system */
  psi= atan2(-yc, -xc);
  cpsi= cos(psi);  spsi= sin(psi);
  gamma= atan(dzds);
  cg=  cos(gamma); sg= sin(gamma);
  /* transform the line  into the helix coordinate system  */
  /* translation */
  lin[0]-= xc;
  lin[1]-= yc;
  lin[2]-= zh;
  /*rotation*/
  x0h=  lin[0]*cpsi +  lin[1]*spsi;
  y0h= -lin[0]*spsi +  lin[1]*cpsi;
  z0h=  lin[2];
  axh=  lin[3]*cpsi +  lin[4]*spsi;
  ayh= -lin[3]*spsi +  lin[4]*cpsi;
  azh=  lin[5];

  /* I did some algebra.. parameterize the line as  (x,y,z)[s] and the
     helix as (x,y,z)[sigma]. Minimize the  distance between them taking
     the restricting case where the line and the z-axis are much much closer
     than 'rho', and solve for sigma.*/

  sigma= ( y0h*cg + z0h*sg + (ayh*cg + azh*sg)*(axh*myrho-axh*x0h-ayh*y0h-azh*z0h)) /
    (1.0 - pow( (ayh*cg + azh*sg),2.0));
  /* at this point the sign of sigma is such that +ve sigma means movement
     along the helix in the positive phi direction (ccw when viewed from 
     downstream). Correct for the   charge of the track so that +ve sigma
     always means moving along the helix in the direction of the particle.*/
  sigma*=sign;

  return sigma;
}
