#include <stdio.h>
#include "rtpc.h"

/////////////////////////////////////////////////////////////
//this routine replaces the old read_paths.c routine and instead
//assigns spatial coordinates for each cell according to the electron
//drift path functions
//////////////////////////////////////////////////////////////

int init_paths()
{
  int pad, side, jj, tt;
 
  printf("gem/init_paths: Initializing RTPC pad locations .....\n");
  generate_pad_locations();

  init_paths_using_fit_eg6();

  // Shift phi by p_path_off:
  for(side=0; side<2; side++)
  {
    for(jj=0; jj<NUM_PADS/2; jj++)
    {
      pad = jj + NUM_PADS/2*side;
      for (tt=0; tt<NAL_SAMP; tt++)
      {
        rawCYL[tt][pad].phi = phrng(rawCYL[tt][pad].phi + p_path_off[side]);
      }
    }
  }

  return 1;

}
