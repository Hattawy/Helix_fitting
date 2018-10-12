#include <stdio.h>
#include <stdlib.h>
#include "map_manager.h"
#include "rtpc.h"
#define NREG 4
//
// March 2014, NAB.  Now reading gains from database.  Divided into 4 groups of
// length 800 because database doesn't like long arrays.  Ordered by increasing
// channel number within each caldb item (previously gain ascii file was ordered
// by row/column).  The four groups are phi-quadrants from beam's perspective:
// (Right/Left/Top/Bottom)
//
//   Channel#     caldb Item
//    0 -  799 == STATUS_RB
//  800 - 1599 == STATUS_RT
// 1600 - 2399 == STATUS_LT
// 2400 - 3199 == STATUS_LB
//
// Internal gem storage of gains array corr_fac is unchanged.
//
int read_scale_db(const int runno)
{
  int dunno,ii;
  char map[200],item[20];
  
  for (ii=0; ii<NUM_PADS; ii++) corr_fac[ii]=1.0;

  // read database, fill corr_fac array:
  printf("Reading RTPC gains from Map/CALDB's RTPC2 for Run# %d ....\n",runno);
  sprintf(map,"%s/Maps/RTPC2.map",getenv("CLAS_PARMS"));
  const char* region[NREG]={"RB","RT","LT","LB"};
  for (ii=0; ii<NREG; ii++)
  {
    sprintf(item,"GAIN_%s",region[ii]);
    if (map_get_float(map,"ENERGY",item,
                      NUM_PADS/NREG,&corr_fac[NUM_PADS/NREG*ii],runno,&dunno))
    {
      fprintf(stderr,"gem/read_scale_db:  Error Reading RTPC2:ENERGY:%s for run %d\n\n",item,runno);
      return 1;
    }
  }

  // count large gain factors:
  int ncrazy=0;
  for (ii=0; ii<NUM_PADS; ii++)
  {
    if (corr_fac[ii] < 0.01)
    {
      if (!dead_pad(ii)) ncrazy++;
      corr_fac[ii]=1.0;
    }
  }
  printf("gem/read_scale_db:  Number of pads (not flagged as bad)  with crazy gain correction factors: %d\n",ncrazy);

  return 0;
}



int read_scale()
     /* read pad by pad gain factors from data file */
{
  char corrfacfilename[200];
  FILE * corrfile;
  int ipad, irow, icol;
  int ncrazy=0;
  typedef char linetype[LINLMAX];
  linetype textline;
  sprintf(corrfacfilename,"pad_gain_factors.dat");

  printf("gem/read_scale:  Reading gain correction file %s\n",corrfacfilename);

  
  if((corrfile = fopen(corrfacfilename,"r"))==NULL)
  {
    printf(">>>>>failed to open gain correction factor file %s\n",corrfacfilename);
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
          return(-1);;
        }
        if(fgets(textline,80,corrfile))
        {
          sscanf(textline, "%g",&corr_fac[ipad]);
          if (corr_fac[ipad] < 0.01) 
          {
            if(!dead_pad(ipad)) ncrazy++;
            corr_fac[ipad]= 1.0;
          }
        }
        else
        {
          printf(">>>>>>> correction factor file too short!\n");
          return(-1);
        }
      }
    }
  printf("gem/read_scale:  Number of pads (not flagged as bad)  with crazy gain correction factors: %d\n",ncrazy);
  return(0);
}
