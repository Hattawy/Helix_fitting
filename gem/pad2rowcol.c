//#include <stdio.h>
#include "rtpc.h"

void pad2rowcol(int pad_num, int*pad_row, int*pad_col)
     /* Get the row and column of pad 'pad_num'
	from knowledge of the pad layout and
	the readout channel order.

     Note: this routine row/col locations to the PADS, not to the 
     ALTRO CHANNELS. Any mistakes in the channel-to-pad correspondence
     should be handled in the data input routine read_BOS.c*/
{
  int rtpc_side, mypad, conn_num, conn_row, conn_col, pad_seq, odd;
  
  /* Synopsis:
     Pads 0-1599 are on the LEFT rtpc; pads 1600-3199 are on the RIGHT rtpc.
     Pad 0 (zero) is near the lower-left-upstream corner of the pad plane.
     Pad 1600 is near the upper-right-upstream corner of the pad plane. 
     
     A 'row' of pads runs along z. It has constant phi. 0-19 on left, 20-39 on right.
     A 'column' of pads runs around the RTPC. It has constant z. 0=upstream, 39=downstream.  */

      rtpc_side= pad_num/1600; /* left=0, right=1 */
      mypad= pad_num % 1600;   /* now 0-1599 */
      
      conn_num= mypad/16;      /* 0-100 connectors on each rtpc*/
      conn_row= conn_num/5;    /* 5 connectors per row */
      conn_col= conn_num % 5;  /* which of the 5 connector columns */
      pad_seq= mypad%16;      /* sequence number within connector: 0-15 */

      /* even/odd connectors use different channel orders */
      odd= conn_num % 2;

      *pad_col=  8*conn_col + seq_col[odd][pad_seq];
      *pad_row=  2*conn_row + seq_row[odd][pad_seq] + 40*rtpc_side;
}
