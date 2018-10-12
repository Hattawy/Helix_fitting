/*Use DVCS solenoid map to return the B-field as a function of radius
  and longitudinal posistion (r,z). Field is returned in kilogauss.

  initial code written by Nate Baillie on 12/03/05
  hcf added r-dependence 12/20/05

prototype:
float bofrz(float,float);
inputs: z position(mm) to find avg field
outputs: the average field (tesla) at the given z is returned */

#include <stdio.h>
#include <math.h>
#include "rtpc.h"

float bofrz(float rr, float zz)
{

  float ss;
  int ii, jj, ll, mm;
  float tt, uu;
  float Bs, Bz, Btot;
  float pointer;
 
  if(zz < -100 || zz > 100) 
    {
    fprintf(stderr,"axial coordinate z is out of range [-100,100] mm: %f\n",zz);
    return(-1);
    }
  if(rr < 0.0 || rr > 60.0)
    {
    fprintf(stderr,"radial coordinate is out of range in bofrz: %f\n",rr);
    return(-1);
    }
 
  ss= rr;

  if(zz<0) zz = -zz;

  pointer = -1000;
  ll = -1000;
  ii=0;
  while((ss > pointer) && (ii<NUM_S1))
    {
    pointer = sVals1[ii];
    if (ss == pointer) ll = ii;
    else if (ss <= pointer) ll = ii-1;
    else ii++;
    }

  pointer = -1000;
  mm = -1000;
  jj=0;
  while((zz > pointer) && (jj<NUM_Z1))
    {
    pointer = zVals1[jj];
    if (zz == pointer) mm = jj;
    else if (zz <= pointer) mm = jj-1;
    else jj++;
    }
  
  tt = (ss - sVals1[ll])/(sVals1[ll+1] - sVals1[ll]);
  uu = (zz - zVals1[mm])/(zVals1[mm+1] - zVals1[mm]);

  Bs = (1-tt)*(1-uu)*int_Bs[ll][mm] + tt*(1-uu)*int_Bs[ll+1][mm] +
    (1-tt)*uu*int_Bs[ll][mm+1] + tt*uu*int_Bs[ll+1][mm+1]; // in Gauss
  Bz = (1-tt)*(1-uu)*int_Bz[ll][mm] + tt*(1-uu)*int_Bz[ll+1][mm] +
    (1-tt)*uu*int_Bz[ll][mm+1] + tt*uu*int_Bz[ll+1][mm+1]; //in Gauss

  Btot = sqrt(pow(Bs, 2)+pow(Bz, 2))/1000; //in kiloGauss
  
  return(Bz/1000.0); /* hcf changed this from Btot to Bz 12/7/05 */
}
