/* This is a subroutine to use the mag field map for the DVCS solenoid to find 
 * the avg field which a track sees at a given z (calculated in another routine 
 * as the midpoint in z between the beginning and end of the track initial code 
 * written by Nate Baillie on 12/03/05

prototype:
inputs: z position(mm) to find avg field
outputs: the average field (kilogauss) at the given z is returned */

//#include <stdio.h>
#include "rtpc.h"

float bofz(float zz)
{

 float ss, xx,yy;
 double bx,by,bz;
 
  if(zz < -100 || zz > 100) 
    {
      errors[3]++;
      return(-1);
    }
  
  ss= 45.0; /* changed to middle of active radial range by hcf 12/7/05 */

  /*
    12/8/06 - hcf attempts to consolidate b-field routines. This routine now 
    calls the more basic and general bfield(x,y,z,&bx,&by,&bz) to obtain the 
    field value.
    Note that we're still assuming radial symmetry wrt the bonus system here, 
    but this is not strictly correct as the solenoid may not be colinear with 
    bonus 
  */

  xx= ss;
  yy= 0.0;
  bfield(xx,yy,zz,&bx,&by,&bz); /* this gives field in Tesla */

  /*bofz is to return kilogauss and sign is expected to be positive*/
  return(-bz*10.0); 
}
