#include "rtpc.h"
//#include <stdio.h>
#include <math.h>

int inside_rtpc(float xx, float yy, float zz)
/* return true if the specified point is within the active rtpc volume */
/* REVISION HISTORY
   9-FEB-07: hcf removed some old commented-out code
  20-MAY-07: hcf reset R_FUZ from 5.01 to 0.01 and created a separate
      routine for use by swim_bonus. swim_bonus needs to be able to swim
      a few mm past the GEM.
   2009 dupre: updated for eg6 run
   2012 nab: added inside_rtpc_cyl
*/

#define Z_FUZ 0.01
#define R_FUZ 0.01
{
  float rr;
  int answer;
  answer=0;
  rr= sqrt(xx*xx + yy*yy);
  if ( (rr >= CAT_RAD-R_FUZ) && 
       (rr <= GEM_RAD+R_FUZ) &&
       (zz >= -(ZEND+Z_FUZ)) &&
       (zz <=  (ZEND+Z_FUZ)) ) answer= 1;
  return(answer);
}

int inside_rtpc_cyl(cylindrical point)
{
  float xx=point.s*cos(point.phi);
  float yy=point.s*sin(point.phi);
  return inside_rtpc(xx,yy,point.z);
}
