#include <stdio.h>
#include "rtpc.h"

int eg6rtpc_read_pedestal()
     /* read pad by pad pedestals data file */
{
  char pedfilename[200];
  FILE * pedfile;
  int ipad, irow, icol;
  int ncrazy=0;
  typedef char linetype[LINLMAX];
  linetype textline;
  sprintf(pedfilename,"eg6rtpc_pedestal.dat");

  printf("gem/read_pedestal:  Reading pedestal file %s\n",pedfilename);

  
  if((pedfile = fopen(pedfilename,"r"))==NULL)
  {
    printf(">>>>>failed to open pedestail file %s\n",pedfilename);
    return -1;
  }
  else
    for(irow=0; irow<80; irow++)
    {
      for(icol=0; icol<40; icol++)
      {
        ipad= rowcol2pad(irow,icol);
        if( (ipad<0) || (ipad>=NUM_PADS) )
        {
          printf(">>>>pad number error!!\n");
          printf(">>>>rowcol2pad(%d,%d) returned %d\n",irow,icol, ipad);
          printf(">>>>aborting\n");
          return(-1);
        }
        if(fgets(textline,80,pedfile))
        {
          sscanf(textline, "%g",&eg6peds[ipad]);
          if (eg6peds[ipad] < 0.0) 
          {
            if(!dead_pad(ipad)) ncrazy++;
            eg6peds[ipad]= 0.0;
          }
        }
        else
        {
          printf(">>>>>>> pedestal file too short!\n");
          return(-1);
        }
      }
    }
  printf("gem/read_pedestal:  Number of pads (not flagged as bad)  with negative pedestals: %d\n",ncrazy);
  return(0);
}
