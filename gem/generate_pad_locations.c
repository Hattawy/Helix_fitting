#include <stdio.h>
#include <math.h>
#include "rtpc.h"
void generate_pad_locations()
     /* Set up the pad locations array
      by creating the locations from knowledge of the pad layout and
     the readout channel order.

     Note: this routine assigns locations to the PADS, not to the 
     ALTRO CHANNELS. Any mistakes in the channel-to-pad correspondence
     should be handled in the data input routine read_BOS.c*/
{
  int pad_num, rtpc_side, mypad, conn_num, conn_row, conn_col, pad_seq, odd;
  int pad_col, pad_row;
  
  /* Synopsis:
     Pads 0-1599 are on the LEFT rtpc; pads 1600-3199 are on the RIGHT rtpc.
     Pad 0 (zero) is near the lower-left-upstream corner of the pad plane.
     Pad 1600 is near the upper-right-upstream corner of the pad plane.
     The X axis (phi=0) extends through the LEFT rtpc.
     The Y axis (phi=PI/2) extends upwards.
  */

  /* Do all the pads... */
  for(pad_num=0; pad_num<3200; pad_num++)
  {
    rtpc_side= pad_num/1600; /* left=0, right=1 */
    mypad= pad_num % 1600;   /* now 0-1599 */
    
    conn_num= mypad/16;      /* 0-100 connectors on each rtpc*/
    conn_row= conn_num/5;    /* 5 connectors per row */
    conn_col= conn_num % 5;  /* which of the 5 connector columns */
    pad_seq= mypad%16;      /* sequence number within connector: 0-15 */

    
    /* even/odd connectors use different channel orders */
    odd= conn_num % 2;
    pad_col=  8*conn_col + seq_col[odd][pad_seq];
    pad_row=  2*conn_row + seq_row[odd][pad_seq];
    
    /* phi_start and zoff give the center coordinates of the corner pads
       one each rtpc module */
    padLocs[pad_num].z= zoff[pad_row%4] + DEL_Z*pad_col;
    padLocs[pad_num].phi= phi_start[rtpc_side] + DEL_PHI*pad_row;
    if(padLocs[pad_num].phi <  0.0) {padLocs[pad_num].phi +=2*PI;}
    if(padLocs[pad_num].phi > 2*PI) {padLocs[pad_num].phi -=2*PI;}

    if (isnan(padLocs[pad_num].z) || isnan(padLocs[pad_num].phi))
      printf("ERROR in generate_pad_locations\n");

/*  printf(">> pad %5d is at %f %f col=%d row=%d\n",
    pad_num,padLocs[pad_num].z, padLocs[pad_num].phi,pad_col,pad_row); */
  }

  //print_pad_locations();

  return;
}

/*
void print_pad_locations()
{
  int row,col,pad;
  for (row=0; row<80; row++)
  {
    for (col=0; col<40; col++)
    {
      pad=rowcol2pad(row,col);
      fprintf(stderr,"row=%2d col=%2d    z=%8.3f phi=%8.3f\n",
          row, col, padLocs[pad].z, padLocs[pad].phi);
    }
  }
}
*/
