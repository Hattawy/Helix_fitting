/*Use DVCS solenoid map to return the B-field as a function of radius
  and longitudinal posistion (r,z). Field is returned in kilogauss.
  initial code written by Nate Baillie on 12/03/05
  hcf added r-dependence 12/20/05
  12/8/06 hcf: incorporate surveyor's offsets between BONUS and DVCS solenoid
prototype:
inputs: z position(mm) to find avg field
outputs: the average field (tesla) at the given z is returned */

#include <stdio.h>
#include <math.h>
#include "rtpc.h"

void bfield(float inxx, float inyy, float inzz, double* bx, double* by, double* bz)
{

  int ii, jj, ll, mm, flip_sign;
  float xx, yy, zz, rr, tt, uu, phi, pointer;
  float Bs, Bz;// /* be aware that Bz != bz !! */
 
  *(bx) = 0;
  *(by) = 0;
  *(bz) = 0;
 
  xx = inxx + XBONUS - XSOL; 
  yy = inyy + YBONUS - YSOL; 
  zz = inzz + ZBONUS - ZSOL;
 
  /* Forgive a 1mm range error in z (just to hush-up the code a little) */
  if (zz < -101.0 || zz > 101.0) return; 
  if (zz < -100.0) zz = -100.0;
  if (zz >  100.0) zz =  100.0;

  rr= sqrt(xx*xx + yy*yy);

  phi = rr>0 ? atan2(yy,xx) : 0.0;

  /* Forgive a 10.1mm out-of-range in r... */
  if (rr>60.0 && rr<70.1) rr=60.0;

  if (rr < 0.0 || rr > 60.0)
  {
    printf("radial coordinate is out of BONUS range in bfield: %f\n",rr);
    errors[2]++;
    return;
  }  

  flip_sign=FALSE;
  if (zz<0) 
  {
    zz = -zz;
    flip_sign=TRUE;
  }

  pointer = -1000;
  ll = -1000;
  ii=0;
  while ((rr > pointer) && (ii<NUM_S1))
  {
    pointer = sVals1[ii];
    if      (rr == pointer) ll = ii;
    else if (rr <= pointer) ll = ii-1;
    else                    ii++;
  }

  pointer = -1000;
  mm = -1000;
  jj=0;
  while ((zz > pointer) && (jj<NUM_Z1))
  {
      pointer = zVals1[jj];
      if      (zz == pointer) mm = jj;
      else if (zz <= pointer) mm = jj-1;
      else                    jj++;
  }
  
  tt = (rr - sVals1[ll])/(sVals1[ll+1] - sVals1[ll]);
  uu = (zz - zVals1[mm])/(zVals1[mm+1] - zVals1[mm]);

  Bs = (1-tt)*(1-uu)*int_Bs[ll][mm] + tt*(1-uu)*int_Bs[ll+1][mm] +
    (1-tt)*uu*int_Bs[ll][mm+1] + tt*uu*int_Bs[ll+1][mm+1]; // in Gauss

  Bz = (1-tt)*(1-uu)*int_Bz[ll][mm] + tt*(1-uu)*int_Bz[ll+1][mm] +
    (1-tt)*uu*int_Bz[ll][mm+1] + tt*uu*int_Bz[ll+1][mm+1]; //in Gauss

  /* CORRECT FOR AN APPARANT SIGN ERROR IN THE FIELD DIRECTION */
  Bz = -Bz;
  Bs = -Bs;
  
  /* Respect radial field sign difference either side of symmetry plane */
  if (flip_sign) Bs = -Bs;

  /* convert to (Bx,By,Bz) in Tesla */
  *(bx)= Bs*cos(phi)/10000.0;
  *(by)= Bs*sin(phi)/10000.0;
  *(bz)= Bz/10000.0;
  
  return; 
}
