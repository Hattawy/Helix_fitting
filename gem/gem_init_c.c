/*      subroutine gem_init_c*/
/*
 *_begin_doc
 *  RCS ID string
 *  $Id: gem_init_c.c,v 1.3.4.4 2013/05/06 15:27:38 claseg6 Exp $
*/
/*
 *  Documentation for subroutine gem_init_c
 *
 *  Purpose:
 *  --------
 *      initialize rtpc data structure (other than geometry)
 *
 *
 *  Calling Sequence:
 *  ----------------
 *
 *  Input Parameters:  (Name - Type - Meaning)
 *  ----------------
 *
 *  Output Parameters:  (Name - Type - Meaning) all in SCS coord system.
 *  -----------------
 
 *  Called from: gem_brun.F 
 *  ------------
 *
 
 *
 *  Notes:  The SCS coordinate system has the nominal target at the origin,
 *  ------  the z-axis along the beamline, the x-axis parallel in the midplane,
 *          y-axis ZxX = pointing up (for Sector 1) parallel to axial wires
 *          pointing towards the HV side (at least for R1).
 *
 *
 * Author: gabriel niculescu, 2005
 *  -------
 *
 * Major revisions:  created - 2005 GN.
 *                   updated for Production run 11-11-05 N. Baillie
 * ----------------

 *
 *_end_doc
 */
/*  RCS information: 
*/

#include <stdio.h>
#include <stdlib.h>
#include "map_manager.h"
#include "rtpc.h"


#define USE(var) static void * use_##var = (void *) &var
static char crcsid[] =
   "$Id: gem_init_c.c,v 1.3.4.4 2013/05/06 15:27:38 claseg6 Exp $";
USE(crcsid);   /* make sure it is not optimized away */
/*
 *   Module information:
 */
  static char crname[] = "GEM_INIT_C";
  static char crauth[] = "Gabriel Niculescu";
/*
 *_begin_inc
 *  include files :
 *  ---------------------
 */
/* BOS common block  uncomment the next line for BOS include file */
/*#include <bos.inc> */
/* _end_inc */

void gem_init_c_()
{
  int ii;
  
//  printf("***************CALLING RTPC INIT ROUTINE GEM_INIT_C*********\n");  
  time (&start);// start = clock();
  //final drift velocity calibration values from special min-I runs
  delta_helium= 0;
  deltaVGL =    0.00;
  deltaVKL =    0;
  deltaVGR =    0.00;
  deltaVKR =    0;
  p_path_off[0] = 0.;
  p_path_off[1] = 0.0;
 
  helium_fraction = 80.0;

  /* Ad-hoc modification of helium from minuit fit */
  printf("gem/gem_init: helium fraction nominally %.2f\n",helium_fraction);
  helium_fraction+= delta_helium;
  printf("gem/gem_init:            ...altered to %.2f\n",helium_fraction);

  /* we read value from database in routine  , should overwrite value below */
  sol_current = 450;//this is the default value  

  nevt_toomanyhits=0;
  frac_nodrift=0;

  qtot_early=0.0, qtot_intime= 0.0, qtot_late= 0.0;
  num_events = 0;
  totTrkd    = 0;
  num_bad_altro_data=0;
  nal_samp = NAL_SAMP;
  nal_keep = NAL_KEEP;
  good_altro_data= FALSE;
  fit_track_to_beamline = 1;
  rwfopt= 1;
  if(rwfopt==0) 
  {
    printf("\n>>>WARNING: FIT TO BEAMLINE WITH RWFOPT=0 MAKES NO SENSE \n");
    printf("-- DUMMY Z-COORDINATE *WILL* BE USED IN FIT!\n");
  }

  //initialize error counters
  for(ii=0; ii<MAX_ERRORS; ii++){errors[ii] = 0;}

  writePadUses = FALSE;
  
  for(ii=0; ii<NUM_PADS; ii++) pad_uses[ii] = 0; 
  
  return;
}

void gem_read_database_(int* irun) {

  //this routine does read the database, but it also initializes some important
  //parameters at the end.

  char map[255];
 // static int firsttime;
 // float caldb_value;
  float VKL;
  float VGL;
  float VKR;
  float VGR;

  /* try to read database value, which is given in mA */
  sprintf(map,"%s/Maps/RUN_CONTROL.map", getenv("CLAS_PARMS") );  
//  map_get_float(map, "currents", "solenoid", 1, &caldb_value, *irun, &firsttime);
//  if (caldb_value != 0) sol_current = caldb_value / -1000.;

//  sol_current == 450;

  fprintf(stdout, "gem/gem_init: Solenoid Current read from database for run %d: %12.3f A\n", *irun, sol_current); 
  
  if(!(sol_current == 450))
    fprintf(stdout, "gem/gem_init: MAJOR WARNING: RTPC is not calibrated for this solenoid current!\n");

  readBFmap();//reads in the solenoid field map.  sol_current must be set
                 //properly in the init_globals routine.

  /* try to read database VGL, which is given in V */
/*  sprintf(map,"%s/Maps/rtpc.map", getenv("CLAS_PARMS") );  
  map_get_float(map, "voltages", "gem_left_rbk", 1, &caldb_value, *irun, &firsttime);
  if (caldb_value != 0) VGL = caldb_value;
  else VGL = -999;*/
  VGL = 2800;

  fprintf(stdout, "gem/gem_init: Voltage Left read from database for run %d: %12.0f V\n", *irun, VGL); 

  /* try to read database VGR, which is given in V */
/*  sprintf(map,"%s/Maps/rtpc.map", getenv("CLAS_PARMS") );  
  map_get_float(map, "voltages", "gem_right_rbk", 1, &caldb_value, *irun, &firsttime);
  if (caldb_value != 0) VGR = caldb_value;
  else VGR = -999;*/
  VGR = 2800;

  fprintf(stdout, "gem/gem_init: Voltage Right read from database for run %d: %12.0f V\n", *irun, VGR);

  /* try to read database VKL, which is given in V */
/*  sprintf(map,"%s/Maps/rtpc.map", getenv("CLAS_PARMS") );  
  map_get_float(map, "voltages", "cathode_left_rbk", 1, &caldb_value, *irun, &firsttime);
  if (caldb_value != 0) VKL = caldb_value;
  else VKL = -999;*/
  VKL = 4300;

  fprintf(stdout, "gem/gem_init: Cathode Voltage Left read from database for run %d: %12.0f V\n", *irun, VKL);

  /* try to read database VKR, which is given in V */
/*  sprintf(map,"%s/Maps/rtpc.map", getenv("CLAS_PARMS") );  
  map_get_float(map, "voltages", "cathode_right_rbk", 1, &caldb_value, *irun, &firsttime);
  if (caldb_value != 0) VKR = caldb_value;
  else VKR = -999;*/
  VKR = 4300;

  fprintf(stdout, "gem/gem_init: Cathode Voltage Right read from database for run %d: %12.0f V\n", *irun, VKR);

  VKL += deltaVKL;
  VGL += deltaVGL;
  VKR += deltaVKR;
  VGR += deltaVGR;

  printf("gem/gem_init: PS voltages after perturbation (VKL, VKR, VGL, VGR) = (%.2f,%.2f,%.2f,%.2f)VOLTS\n",
	 VKL, VKR, VGL, VGR);

  calcVolt(VGL, VGR, VKL, VKR);  
  runNum = *irun;
  init_paths();  //this routine initializes our time to space coordinate conversion
  altro_offset= 16.0;
  altro_thresh=5.0;
  init_CS();


  // READ FROM CALIBRATION DATABASE (March 2014, NAB): 
  if (read_bad_pads_db(*irun)) exit(0);
  if (read_scale_db(*irun))    exit(0);

  // READ FROM ASCII FILES:
  //if (read_bad_pads()) exit(0);
  //if (read_scale())    exit(0);

  // Get pedestals:
  //if(eg6rtpc_read_pedestal()) {printf(">>>>> error return from read_pedestal\n"); exit(0);}
 
}
