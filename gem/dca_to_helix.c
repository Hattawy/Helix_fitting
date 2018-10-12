#include "rtpc.h"
#include <stdio.h>
#include <math.h>

float dca_to_helix(float rho, float x_c, float y_c, float dzds, float x_s ,
                   float y_s, float z_s, hitlist *apoint, cylindrical* close, 
                   float *sigma_d) 
/* Given 1. A helix specified by
      rho -- radius (-ve for -ve track)
      x_c  = x of center
      y_c  = y of center
      dzds = tan(angle out of xy plane)
      x_s  = x of starting point of helix
      y_s  = y  "
      z_s  = z  "

 2. An arbitrary space point specified in hitlist coordinates by
      (*apoint).r
      (*apoint).phi
      (*apoint).z

Find the point on the helix closest to the point and return it as
      close.s
      close.phi
      close.z
         and
      float *sigma = path length along helix from (x_s,y_s,z_s) to that point.

Return value of the function is the distance from the helix to the point at
the point of closest approach.  */
{

  float gamma, cg, sg, x_p, y_p,z_p, phi_p, phi_0, del, sigma, x_d, y_d, z_d;
  float phi_ph, z_ph, x_dh, y_dh, z_dh, delphi, delz,  radius, num_turns;
  int sign_q;

  if(GEMDEBUG) {printf("dca_to_helix given (rho,x_c,y_c)= (%5.2f,%5.1f,%5.1f), dzds=%5.1f, (x_s,y_s,z_s)= (%5.1f,%5.1f,%5.1f)\n",
		       rho,x_c,y_c, dzds, x_s,y_s,z_s);}

  gamma= atan(dzds);  /*dz is dz. ds is projection of path length into xy plane.*/
  sg= sin(gamma);
  cg= cos(gamma);
  radius= fabs(rho);
  sign_q= rho/radius;

  /* some things are just easier in cartesian coordinates */
  /* the arbitrary point... */
  x_p= (*apoint).r*cos((*apoint).phi);
  y_p= (*apoint).r*sin((*apoint).phi);
  z_p= (*apoint).z;

  /* Relative to helix axis, find xy angle to starting point of helix */
  phi_0= atan2(y_s-y_c, x_s-x_c);

  /* Similarly for the space point */
  phi_p= atan2(y_p-y_c, x_p-x_c);

  /* The angle of the point within the helix's own coordinate system */
  /* Note that there is an arbitrary offset of 2*PI*N where N is undetermined */
  phi_ph= fabs(delphrng(phi_p - phi_0));

  /* Z of the arbitrary point in the helix system...*/
  z_ph= (*apoint).z - z_s;

  /* Now we can find the helical path length to the closest point... derived
     by some simple calculus. Aftwards we attempt to correct for the
     arbitrary 2*PI*N in the solution by asking for a reasonable Z */

  
  sigma= radius*phi_ph*cg + z_ph*sg;

  /* Given sigma, get the coordinates of that point on the helix. */
  /* In the helix coordinate system: */
  x_dh= radius*cos( sigma*cg/rho );
  y_dh= radius*sin( sigma*cg/rho );
  z_dh= sigma*sg;

  /* Transform back to the lab coordinate system */
  x_d = x_c + x_dh*cos(-phi_0) + y_dh*sin(-phi_0);
  y_d = y_c - x_dh*sin(-phi_0) + y_dh*cos(-phi_0);
  z_d = z_dh + z_s;

  /* And, finally, into lab cylindrical coordinates */
  (*close).s  = sqrt(x_d*x_d + y_d*y_d);
  (*close).phi= phrng(  atan2(y_d, x_d) );
  (*close).z  = z_d;

  /* Calculate distance between space point and nearest point on helix: */
  delz= (z_p-z_d);
  delphi= delz/(radius*dzds);
  /* if |delphi| is > PI then we probably need to add or subtract 2PI. */
  num_turns= (float)rint(delphi/(TWOPI));
  /* If we took off or added a whole turn (or more) recompute the coordinates */
  /* but do not add or remove turns if the delz is less than our z resolution */
  if(  (fabs(num_turns)>0.1) && (fabs(delz)>7.0) )
  {
    phi_ph+= 2.0*PI*num_turns;
    sigma= radius*phi_ph/cg;
    /* Given sigma, get the coordinates of that point on the helix. */
    /* In the helix coordinate system: */
    x_dh= radius*cos( sigma*cg/rho );
    y_dh= radius*sin( sigma*cg/rho );
    z_dh= sigma*sg;

    /* Transform back to the lab coordinate system */
    x_d = x_c + x_dh*cos(-phi_0) + y_dh*sin(-phi_0);
    y_d = y_c - x_dh*sin(-phi_0) + y_dh*cos(-phi_0);
    z_d = z_dh + z_s;

    /* And, finally, get the point on the helix into lab cylindrical coordinates */
    (*close).s  = sqrt(x_d*x_d + y_d*y_d);
    (*close).phi= phrng(atan2(y_d, x_d));
    (*close).z  = z_d;

  }

  del= sqrt( (x_p-x_d)*(x_p-x_d) + (y_p-y_d)*(y_p-y_d) + (z_p-z_d)*(z_p-z_d) );

  /* Store the residuals to each point */
  (*apoint).res_r= (*apoint).r - (*close).s;
  (*apoint).res_phi= delphrng( (*apoint).phi - (*close).phi );
  (*apoint).res_z= (*apoint).z - (*close).z;

   if(GEMDEBUG) {printf("nearest point on helix to (%5.1f,%5.1f,%5.1f) is (%5.1f,%5.1f,%5.1f), del=%5.2f sigma=%5.1f\n", 
 		       x_p,y_p,z_p, x_d,y_d,z_d, del, sigma);}        

  *sigma_d= sigma;
  return del;
}
