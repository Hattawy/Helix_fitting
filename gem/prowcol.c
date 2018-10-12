#include "rtpc.h"
//#include <stdio.h>
 
int prowcol(float zz, float phi, float *prow, float *pcol)
     /* given (zz,phi) coordinates on the pad plane, return
	prow= row number and
	pcol= column number of the pad covering that coordinate, and
	function value= the actual pad number of that pad (0-3199)
	   return -1 if electron did not drift to a pad.

        phi should be in radians, with phi=0 defined as being
	in the middle of the LEFT RTPC.

	Pad row 0 (zero) is at the lower-left of the rtpc when
	viewed from upstream, 39 is upper-left, 40 is upper-right,
	and 79 is lower right.  */
{
  float myphi, zstart;
  int mycol, myrow, rtpc_side, conn_row, conn_side, conn_col, conn_num, even_odd, padnum ;
  int rows_over;

  /*rtpc coordinates have phi=0 along x axis, but pad numbering
    starts at bottom-upstream of left module.  */

   myphi= phi;
   rows_over= 0;
  for(;myphi<0.0; myphi+=2.0*PI);
  for(;myphi>2.0*PI; myphi-=2.0*PI);
  
  if( (myphi>0.0) && (myphi<PI/2.0) )
  {
    rtpc_side=0; /* LEFT RTPC module */
    rows_over= (myphi-0.0)/DEL_PHI;
    myrow= 20 + rows_over;
  }
  else if ( (myphi>3.0*PI/2.0) )
  {
    rtpc_side=0; /* LEFT RTPC module */
    rows_over= (myphi-(phi_start[rtpc_side] - (DEL_PHI/2.0)))/DEL_PHI;
    myrow=  0 + rows_over;
  }
  else
  {
    rtpc_side=1; /* RIGHT RTPC module */
    rows_over= (myphi-(phi_start[rtpc_side] - (DEL_PHI/2.0)))/DEL_PHI;
    myrow= 40 + rows_over;
  }

  if( (myrow<0) || (myrow>79) || (rows_over<0) ) return(-1);

  *prow= myrow; /* return this */
  /* Get the local (per side) index of this row */
  myrow= myrow % 40;
  /* Now, myrow should be in the range 0-39 */

  /* The z of the edge of the first pad on this row... */
  zstart= zoff[myrow%4] - DEL_Z/2.0;

  /* Which of the 40 columns (along phi) is this pad in?*/
  mycol= (zz-zstart)/DEL_Z;

  /* There are 20 rows of connectors, each one picking up two rows of pads */
  conn_row= myrow/2; 

  /* Within the connector there are two sides */
  conn_side= myrow %2;

  /* ...and each side sees 8 pads */
  conn_col= mycol % 8;

  /* The connectors are organized as 5 per row */
  conn_num= 5*conn_row + mycol/8;

  /* The channels within each connector are ordered in a way to make the pcb
     layout convenient, so the order appears to be messy. The order is different
     for even and odd-numberd connectors, as these correspond to ALTRO chips on
     opposite sides of the FEC pcb. */

  even_odd= conn_num %2;

  /* Finally we can calculate the pad number = ALTRO channel number */
  if( (myrow>=0) && (myrow<40) && (zz>=zstart) && (mycol<40) )
    padnum= 16*conn_num + pcb_map[even_odd][conn_side][conn_col] + 1600*rtpc_side;
  else
    padnum=-1;

  *pcol= mycol;

  return(padnum);
}
