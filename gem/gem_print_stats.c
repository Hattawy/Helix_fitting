#include "rtpc.h"
#include <stdio.h>
#include <time.h>
void gem_print_stats()
{
  // Moved all this to get it out of gem_evnt. (nab 2012)

  if(num_events == 1)
  {	    
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf("Begin first RTPC reconstruction. Current time: %s",asctime(timeinfo));
  }
  else if (!(num_events % 1000)) 
  {
    time (&end);
    elapsed = difftime (end,start);
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf("%d TPC trigs analyzed in %.2f minutes. Current time: %s",
        num_events,elapsed/60.0,asctime(timeinfo));
    printf("Avg = %.2f min / 10000 events\n\n",10000*elapsed/(float)num_events/60.0);
  }
}
