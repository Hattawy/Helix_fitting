/* this routine is called once at the beginning of main to initialize the 
 * rotation angles between the various coordinate systems, all measurements 
 * are taken wrt the CLAS CS.  more details can be found in the documentation 
 * of the translation routines transXYZ and transMOM
 *
 * Below the results of the Hall B survey should be implemented
 * 
 * Author: Nate Baillie
 * Date Created: 12-14-06 */
#include <stdio.h>
#include <math.h>
#include "rtpc.h"

void init_CS()
{
  double bmx0_here, bmdxdz_here, bmy0_here, bmdydz_here;
  GetBLPara(runNum, 1, &bmx0_here, &bmy0_here, &bmdxdz_here, &bmdydz_here);
  bmx0 =   bmx0_here;
  bmy0 =   bmy0_here;
  bmdxdz = bmdxdz_here;
  bmdydz = bmdydz_here;
  //I want the x and y intercepts at Z= ZBONUS ~ -580*/
  bmx0= bmx0 + bmdxdz*(ZBONUS+640.);
  bmy0= bmy0 + bmdydz*(ZBONUS+640.);
  
  if(runNum <=10)
  {
    printf("this appears to be simulated data, using ideal beamline\n");
    bmx0 =   0.0;
    bmy0 =   0.0;
    bmdxdz = 0.0;
    bmdydz = 0.0;
  }
  printf("*****beam position parameters used for run %d=\n",runNum);
  printf("*****bmx0=%.5lf bmy0=%.5lf bmdxdz=%.5lf bmdydz=%.5lf \n",
	 bmx0, bmy0, bmdxdz,bmdydz );
  bmalpha = 0.0;
  bmbeta  = (float)atan(bmdydz);
  bmgamma = (float)atan(bmdxdz);
  bmxoff  = (float)(bmx0 + bmdxdz*640.0);
  bmyoff  = (float)(bmy0 + bmdydz*640.0);
  bmzoff  = 0.0;

  boalpha = 0.0;
  bobeta  = 0.0;
  bogamma = 0.0;
  boxoff  = 0.0;
  boyoff  = 0.0;
  bozoff  = -630.0;

  soalpha = 0.0;
  sobeta  = 0.0;
  sogamma = 0.0;
  soxoff  = 0.0;
  soyoff  = 0.0;
  sozoff  = -630.0;
  return;
}
