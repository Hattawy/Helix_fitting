//#include <stdio.h>
#include "rtpc.h"
 
int rowcol2pad(int prow, int pcol)
     /* Given the row and column of a pad, figure out its pad number.
	Return the pad number as the function value.

	This is a complementary function to prowcol.c. 

	It is intended to be used to get the pad number corresponding
	to a pad's location on a row vs. column scatter plot */
{
  int myrow, mycol;
  int even_odd, conn_num, conn_row, conn_col, conn_side, rtpc_side, padnum;

  myrow= prow;
  mycol= pcol;

  /* There are two rtpc modules... */
  rtpc_side= myrow/40;
  myrow= myrow % 40;

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
  if( (myrow>=0) && (myrow<40) )
    padnum= 16*conn_num + pcb_map[even_odd][conn_side][conn_col] + 1600*rtpc_side;
  else
    padnum=-1;

  return(padnum);
}

