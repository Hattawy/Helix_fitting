#include <stdio.h>
#include <stdlib.h>
#include "rtpc.h"

/*this routine will read in the solenoid field map given to us by saclay.  it fills
  in a couple arrays int_Bs, and int_Bz and passes them back to bofz so we can 
  determine the average B field at a given z*/

void readBFmap()
{
  FILE *magFile;
  
  int ii, jj, kk;
  float temp_s[NUM_VAL1];
  float temp_z[NUM_VAL1];; //the array of the s and z values from the field map
  float temp_Bs[NUM_VAL1];
  float temp_Bz[NUM_VAL1];
  float trash;
  float CurrentRatio;
  // CurrentRatio is: expect B at origin/the B in data map at origin
  // In Bonus case, CurrentRatio=2.0/4.6884=0.426585
  // changed for 5/11/05 simulation to 2.5/4.6884 = 0.533230953
  // for CLASrun we ran with lots of different B fields - the ratio 
  // for currents scales like (current setting)/ 534A
  // for 400A -> .7490636704

  CurrentRatio = sol_current/534.0;
  
  char filename[1000];
  sprintf(filename,"%s/eg6rtpcBFieldMap.txt",getenv("CLAS_PARMS"));

  if((magFile = fopen(filename,"r"))==NULL)
//  if((magFile = fopen("BFieldMap.txt","r"))==NULL)
  {
    printf("gem/readBFmap: failed to open magnetic field data file\n");
    exit(1);
  }
  else 
  {
    fgets(textline, LINLMAX, magFile); //eat the first line
    printf("gem/readBFmap: read magnetic field file\n");
  }
  for (ii=0; ii<NUM_VAL1; ii++)
  {
    fscanf(magFile,"%f %f %f %f %f %f",&trash,&temp_s[ii],&temp_z[ii], &trash, &temp_Bs[ii], &temp_Bz[ii]);
    /* we just read the s, z and Bs, Bz values*/
    //printf("read s = %f z = %f Bs = %f Bz = %f\n", temp_s[ii], temp_z[ii], temp_Bs[ii], temp_Bz[ii]);
  }
  jj = 0;
  for (ii=0; ii<NUM_VAL1; ii+=NUM_Z1)
  {
    sVals1[jj]=10*temp_s[ii]; //convert to mm
/*     printf("sVals1[%d] = %f \n", jj, sVals1[jj]); */
    jj++;
  }
  for (ii=0; ii<NUM_Z1; ii++)
  {
    zVals1[ii] = 10*temp_z[ii] + 790; //put origin at z = 0 and convert to mm
/*     printf("zVals1[%d] = %f \n", ii, zVals1[ii]); */
  }
  kk=0;
  for (ii = 0; ii < NUM_S1; ii++)
  {
    for (jj = 0; jj < NUM_Z1; jj++)
    {
      int_Bs[ii][jj] = temp_Bs[kk] * CurrentRatio;
      int_Bz[ii][jj] = temp_Bz[kk] * CurrentRatio; 
/*    printf("s and z %f %f gives int_Bs = %f and int_Bz = %f \n", sVals1[ii], zVals1[jj], int_Bs[ii][jj], int_Bz[ii][jj]); */
      kk++;
    }
  }
}
