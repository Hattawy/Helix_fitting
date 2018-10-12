#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "map_manager.h"
#include "rtpc.h"
#define NREG 4
//
// March 2014, NAB.  Now reading pad status flags from database, and badPads array
// is now just a copy of those 3200 statuses (0=good, 1=bad), while previously it
// contained bad channel numbers.  Divided into 4 caldb items of length 800 because
// database doesn't like long arrays.  Ordered by increasing channel number within each
// caldb item.   The four groups are phi-quadrants from beam's perspective:
// (Right/Left/Top/Bottom)
//
//   Channel #    caldb Item
//    0 -  799 == STATUS_RB
//  800 - 1599 == STATUS_RT
// 1600 - 2399 == STATUS_LT
// 2400 - 3199 == STATUS_LB
//
// dead_pad() and read_bad_pads() are updated for new format of badPads array
//
int read_bad_pads_db(const int runno)
{
  int dunno,ii;
  char map[200],item[20];
  
  for (ii=0; ii<NUM_PADS; ii++) badPads[ii]=0;

  // read database, fill badPads array:
  printf("Reading RTPC pad statuses from Map/CALDB's RTPC2 for Run# %d ....\n",runno);
  sprintf(map,"%s/Maps/RTPC2.map",getenv("CLAS_PARMS"));
  const char* region[NREG]={"RB","RT","LT","LB"};
  for (ii=0; ii<NREG; ii++)
  {
    sprintf(item,"STATUS_%s",region[ii]);
    if (map_get_int(map,"FLAGS",item,
                    NUM_PADS/NREG,&badPads[NUM_PADS/NREG*ii],runno,&dunno))
    {
      fprintf(stderr,"gem/read_bad_pads_db:  Error Reading RTPC2:FLAGS:%s for run %d\n\n",item,runno);
      return 1;
    }
  }

  // count bad pads:
  int nbad=0;
  for (ii=0; ii<NUM_PADS; ii++)
    if (badPads[ii]==1) nbad++;
  printf("gem/read_bad_pads_db: Number of pads flagged as bad: %d = %.1f%%\n",nbad,100*(float)nbad/NAL_CHAN);
  
  return 0;
}


// Read bad pads from ascii file:
void read_bad_pads()
{
  int done, ii, holder;
  char comment[200], textline[200];
  FILE *badPadfile;

  // default all pads are good:
  for (ii=0; ii<NAL_CHAN; ii++) badPads[ii]=0;

//  float prow, pcol;
#define NODATA -1000
  if((badPadfile = fopen("badPads.dat","r"))==NULL)
  {
    printf("gem/read_bad_pads:  failed to open the hot pad file\n");
    exit(1);
  }
  else
  {
    done = FALSE;
    ii=0;
    while(!done)
    {
      if(ii>NAL_CHAN) // hotPads is dimensioned NAL_CHAN
      {
        printf(">>>ERROR: List overflow in read_bad_pads!! Aborting.\n");
        exit(0);
      }
      if(fgets(textline,200,badPadfile) == NULL)
      {
        done = TRUE;
        num_bad_pads = ii;
      }
      else
      {
        sprintf(comment,"no comment\n");
        sscanf(textline, "%d %100c", &holder,comment);
        if (holder == NODATA) 
        {
          done = TRUE;
          num_bad_pads = ii;
        }
        else 
        {
          comment[strlen(comment)-1]= 0;  // delete newline character
//          printf("GEM: suppressing signals from bad channel no. %d {%s}\n",holder,comment);

	  if (holder<0 || holder>=NAL_CHAN)
	  {
	    fprintf(stderr,"gem/read_bad_pads  ERROR, invalid pad#:  %d\n",holder);
	    exit(1);
	  }

// March 2014, changed format of badPads array for caldb, and 
// updated this routine for backgward compatability (so ascii file still works):
//        badPads[ii] = holder; // OLD
	  badPads[holder] = 1; // NEW


	  //          if(!(prowcol(padLocs[holder].z,padLocs[holder].phi, &prow, &pcol)>=0))
//            printf("Bad pad number found in read_bad_pads: %d\n",holder);
        }
        ii++;
      }
    }
    printf("gem/read_bad_pads: Number of pads flagged as bad: %d = %.1f%%\n",num_bad_pads,100*(float)num_bad_pads/NAL_CHAN);
  }
} 
