#include "rtpc.h"

// March 2014, After changing to use caldb (read_bad_pads_db),
// badPads array now stores 3200 statuses (1=bad, 0=good).  So
// stop searching array, just access one element.
int dead_pad(int pad_num)
{
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // turn back on a few boards!!!!!!!!!!!!!!!!!
  const int brd = pad2board(pad_num);
  if (brd==19 || brd==38) return FALSE;
  if (brd==73 || brd==74) return FALSE;
  if (brd==124 || brd==134) return FALSE;

  // return TRUE if invalid or bad pad_num, else FALSE for good pad
  if (pad_num < 0 || pad_num >= 3200) return TRUE;
  if (badPads[pad_num]==1) return TRUE;
  else                     return FALSE;
}




//
// March 2014, changed badPads array format, this function is dead:
//int dead_pad(int pad_num)
//     /* Returns TRUE if pad_num is in the hot/dead pads list 
//        or is an invalid pad number. */
//{
//  int ii;
//  if ((pad_num < 0) || (pad_num > 3200) ) return (TRUE);
//  for (ii=0; ii<num_bad_pads; ii++)
//  {
//    if(pad_num == badPads[ii]) return (TRUE);
//  }
//  return(FALSE);
//}
//

