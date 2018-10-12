#include "rtpc.h"

void z2beamxy(float z, float*xx, float*yy)
     /* given a Z position in BONUS coordinates, return the (x,y)
	position in BONUS coordintes of the beam at that z*/
{
  *xx = bmdxdz * z + bmx0 - XBONUS;
  *yy = bmdydz * z + bmy0 - YBONUS;

  return;
}
