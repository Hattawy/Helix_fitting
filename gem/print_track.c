#include "rtpc.h"
#include <stdio.h>
#include <math.h>
 
void print_track(int ntracks)
{
  float p_perp, track_theta_deg, track_phi_deg;

  p_perp= sqrtf(pow(trackv[ntracks].p_x,2)+pow(trackv[ntracks].p_y,2));
  track_theta_deg = trackv[ntracks].theta*180.0/PI;
  track_phi_deg= track_phi*180.0/PI;

  fprintf(stderr,"------------------------------------------------------\n");
  fprintf(stderr,"Fit Results for Track %d (Chain %d) eid=%d\n",ntracks, track_chain[ntracks],clas_global_eid);    
/*-----------------------------------------------------------------*/
  fprintf(stderr,"| Ptot: %7.1f  ",trackv[ntracks].sign * trackv[ntracks].p_tot);
  fprintf(stderr,"| THTA: %7.1f  ",track_theta_deg);
  fprintf(stderr,"|  dca: %7.1f  ",dca[ntracks]);
  fprintf(stderr,"|  vtl: %7.1f  ",trackv[ntracks].vtl);

  fprintf(stderr,"\n");
  /*-----------------------------------------------------------------*/
  fprintf(stderr,"|   Px: %7.1f  ",trackv[ntracks].p_x);
  fprintf(stderr,"|  PHI: %7.1f  ",track_phi_deg);
  fprintf(stderr,"| dqdx: %7.0f  ",trackv[ntracks].dqdx);
  fprintf(stderr,"| npts: %7.0f  ",trackv[ntracks].npts);
    

  fprintf(stderr,"\n");
  /*-----------------------------------------------------------------*/
  fprintf(stderr,"|   Py: %7.1f  ",trackv[ntracks].p_y);
  fprintf(stderr,"| thta: %7.3f  ",trackv[ntracks].theta);
  fprintf(stderr,"| edis: %7.1f  ",trackv[ntracks].edist);
  fprintf(stderr,"|    Q: %7.0f  ",trackv[ntracks].visq);


  fprintf(stderr,"\n");
  /*-----------------------------------------------------------------*/
  fprintf(stderr,"|   Pz: %7.1f  ",trackv[ntracks].p_z);
  fprintf(stderr,"|  phi: %7.3f  ",trackv[ntracks].phi);
  fprintf(stderr,"| sdis: %7.1f  ",trackv[ntracks].sdist);


  fprintf(stderr,"\n");
  /*-----------------------------------------------------------------*/
  fprintf(stderr,"|  Pxy: %7.1f  ",p_perp);
  fprintf(stderr,"|  Z_v: %7.1f  ",trackv[ntracks].zdca);
  fprintf(stderr,"| c2/d: %7.2f  ",trackv[ntracks].chi2);
  fprintf(stderr,"| lfit: %7.2f  ",trackv[ntracks].fitqual);


  fprintf(stderr,"\n");
  /*-----------------------------------------------------------------*/
  fprintf(stderr,"| Pcor: %7.1f  ",trackv[ntracks].sign * trackv[ntracks].p_corr);
  fprintf(stderr,"|   R0: %7.1f  ",trackv[ntracks].r_0);

  fprintf(stderr,"\n");
  
  return;
}
