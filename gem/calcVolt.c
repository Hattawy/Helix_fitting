#include <stdio.h>
#include "rtpc.h"

////////////////////////////////////////////////////////////////
//This routine calculates the actual potential difference for 
//each chamber half given the voltage settings on the cathode
//and the first GEM.  This routine puts the spreadsheet we usually
//use for this purpose into c form. 
//2-28-07 edited so that we calculate the drift voltages AND the
//        voltage on the INNER side of the FIRST gem.
////////////////////////////////////////////////////////////////

void calcVolt(float VGL, float VGR, float VKL, float VKR)
{
  float R1, R2, R3, R4, R5, I1, I3;

  printf("gem/calcVolt: nouvelles resistances actives\n");
  //start with the left
  R1  = 2.00e+06; //VK limiter in spread sheet
  R2  = 0.99e+06; //VG limiter
  R3  = 7.32e+07; //sum of field cage resistors 
  R4  = 1.42e+07; //sum of gem resistors
  R5  = 1.27e+07; //same than 4 excluding the first resistor on gem 1
  I1  = (VKL-VGL+(VGL*R2)/(R4+R2))/(R2-(R2*R2)/(R4+R2)+R3+R1);
  I3  = (I1*R2+VGL)/(R4+R2);
  Lvdrift = I1*R3;
  Lvgem   = I3*R5;

  //now the right
  R1  = 2.00e+06;
  R2  = 1.00e+06;
  R3  = 7.35e+07;
  R4  = 1.43e+07;
  R5  = 1.28e+07;
  I1  = (VKR-VGR+(VGR*R2)/(R4+R2))/(R2-(R2*R2)/(R4+R2)+R3+R1);
  I3  = (I1*R2+VGR)/(R4+R2);
  Rvdrift = I1*R3;
  Rvgem   = I3*R5;

  printf("gem/calcVolt:  (Lvdrift, Rvdrift, Lvgem, Rvgem) = (%.0f,%.0f,%.0f,%.0f)VOLTS\n",
	 Lvdrift,Rvdrift,Lvgem,Rvgem);
}
