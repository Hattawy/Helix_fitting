#include <stdio.h>
#include <math.h>
#include "rtpc.h"
// nab modified 2012

// called at end of events -- in recsis/rendrn.F and renjob.F:
int gem_erun_()
{
  int ii;

  // scaling these by the total number of tracks:
  float percenterr[MAX_ERRORS];
  for(ii=0;ii<MAX_ERRORS;ii++) 
    percenterr[ii] = 100*(float)errors[ii]/(float)totTrkd;

  // scaling these differently:
  qtot_early/= (float)num_events*(TPC_TGEM);
  qtot_intime/= (float)num_events*(1+TPC_TCAT-TPC_TGEM);
  qtot_late/= (float)num_events*(nal_samp-TPC_TCAT);

  printf("\n----------------------------------- GEM INFO -------------------------------\n");
  gem_print_stats();
  printf("**** %d events skipped because of not good altro data\n",num_bad_altro_data);
  
  printf(" Q_EARLY  : %.3f per event*time_bin\n",qtot_early);
  printf(" Q_INTIME : %.3f per event*time_bin\n",qtot_intime);
  printf(" Q_LATE   : %.3f per event*time_bin\n",qtot_late);

  printf("****RTPC reconstruction report****\n");
  printf("SUBROUTINE        INCIDENT                                    OCCURED     %% OF TRACKS\n");
  printf("rwfthc()          called w/ > 300 points                     %8d       %.2f\n",errors[0],percenterr[0]);
  printf("rwfthc()          called w/ < 3 points                       %8d       %.2f\n",errors[1],percenterr[1]);
  printf("bfield()          radial coordinate is out of range          %8d       %.2f\n",errors[2],percenterr[2]);
  printf("bofz()            axial coordinate is out of range           %8d       %.2f\n",errors[3],percenterr[3]);
  printf("chain_link()      given too many chains                      %8d       %.2f\n",errors[4],percenterr[4]);
  printf("helix_fit()       called with wrong chain number             %8d       %.2f\n",errors[5],percenterr[5]);
//  printf("read_BOS()        more hits than space allocated in hitlist  %8d       %.2f\n",errors[6],percenterr[6]);
  printf("read_BOS()        channel from bos bank out of range         %8d       %.2f\n",errors[7],percenterr[7]);
  printf("read_BOS()        altro time value not allowed               %8d       %.2f\n",errors[8],percenterr[8]);
  printf("read_BOS()        altro time value out of range              %8d       %.2f\n",errors[9],percenterr[9]);
  printf("read_BOS()        more extras than space allocated           %8d       %.2f\n",errors[10],percenterr[10]);
  printf("scale_hits()      finds a bad pad number                     %8d       %.2f\n",errors[11],percenterr[11]);
  printf("sum_hit_charges() finds a pad number disagreement            %8d       %.2f\n",errors[12],percenterr[12]);
  printf("track_charge()    track never intersects RTPC                %8d       %.2f\n",errors[13],percenterr[13]);
  printf("track_charge()    track starts inside RTPC never hits a wall %8d       %.2f\n",errors[14],percenterr[14]);
  printf("%% of events with too many hits (truncated)                            %.2f\n",100*(float)nevt_toomanyhits/num_events);
  printf("%% of tracks using bad TDC                                             %.2f\n",100*(float)frac_nodrift/totTrkd);
  printf("--------------------------------END GEM INFO -------------------------------\n\n");

  
/*
  if(writePadUses)
  {
    FILE *puf; //pad use file
    char pufname[200];
    sprintf(pufname,"/u/group/bonus/nate/pad_uses/pad_uses%d.dat",runNum);
    
    if((puf = fopen(pufname,"w"))==NULL)
    {
      printf("failed to open %s\n",pufname);
    }
    else
    {
      for(ii=0; ii<NUM_PADS; ii++)
      {
        fprintf(puf,"%d \n",pad_uses[ii]);
      }
    }
    fclose(puf);
  }
*/

  return 0;
}
