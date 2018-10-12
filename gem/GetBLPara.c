#include <stdio.h>
#include <stdlib.h>
#include "map_manager.h"

// nab modified 2012

void GetBLPara(int run, int opt, double *x0, double *y0, double *a, double *b )
{//input: 
// run: run number
// opt:(option) 0: Get Raw Value; "none zero": Get Fitted Value
// output
// beam position parameters: x0,y0,a,b, where X=x0+a*(z+590.),  Y=y0+b*(z+590.) 

/*
   typedef struct {
   int Run;
   double X0;
   double Y0;
   double A;
   double B;
   int Opt;
   }BeamLine;
   BeamLine thisRunBL;
   thisRunBL.Run = -1;
   thisRunBL.X0 = 0.;
   thisRunBL.Y0 = 0.;
   thisRunBL.A = 0.;
   thisRunBL.B = 0.;
   thisRunBL.Opt = 1;

   run = 61000;
*/

  int ii,firsttime;
  float xyz[3];
  char map[255];

  *x0 = 0.;
  *y0 = 0.;
  *a = 0.;
  *b = 0.;
  opt = 1;


/*
// OLD BEAM/RTPC OFFSETS

//const float ebeam[5]   = {     5.7,     1.2,   5.7,    6.0,    1.2     };
  const int runranges[6] = {61001 , 61225 , 61483 , 61510 , 61932 , 61972};
  const float xoffset[5] = {      0,     1.55,   2.37,   2.37,    0      };
  const float yoffset[5] = {      0,     0.29,  -0.4,   -0.4,     0      };
//const float ebeam[5]   = {     5.7,     1.2,   5.7,    6.0,    1.2     };

  if (run < runranges[0] || run >= runranges[5])
  {
    fprintf(stderr,"gem/GetBLPara don't know about this run#: %d\n",run);
    return;
  }

  for (ii=0; ii<5; ii++)
  {
    if (run>=runranges[ii] && run<runranges[ii+1])
    {
      *x0 = xoffset[ii];
      *y0 = yoffset[ii];
      break;
    }
  }
*/

// NEW BEAM/RTPC OFFSET
*x0=0.28;
*y0=0.97;

/*
  sprintf(map,"%s/Maps/GEOMETRY.map",getenv("CLAS_PARMS"));
  map_get_float(map,"beam","position",3,xyz,run,&firsttime);
  fprintf(stderr,
      "gem/GetBLPara.c got beam position from CLAS_PARMS/Maps/GEOMETRY:  %.1f,%.1f,%.1f\n",
      xyz[0],xyz[1],xyz[2]);

  *x0=xyz[0];
  *y0=xyz[1];
*/

  return;
}


