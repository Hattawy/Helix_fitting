#include <stdio.h>
#include <ntypes.h>
#include <bostypes.h>
#include <bosddl.h>     //clasTPC_t,tpc_t was defined here
#include "bosio.h"
#include "bosfun.h"
#include "rtpc.h"
//
// nab modified 2012
//

#define READBOS_DEBUG 0
//global variabls
int  clas_event_type=0;
int  clas_trig_type=0;

int fill_hits(clasTPC_t *pTPC, int bankindex);

int read_BOS()
{
  int i=0,j=0,readflag;
  //clasTPC_t *TPC1=NULL;
  //clasTPC_t *TPC2=NULL;
  clasHEAD_t *HEAD=NULL;
  clasTGBI_t *TGBI=NULL;
  
  //get the eventid from the head bank 
  HEAD = (clasHEAD_t*) getBank(&bcs_, "HEAD");
  if(HEAD)
  {
    clas_global_eid =(int)HEAD->head[0].nevent;
    clas_event_type=(int)HEAD->head[0].type; //if less than 0, it is simulation event
  }
  else
  {//if No HEAD bank, no need to read this event
    clas_global_eid=-1;
    return -1;
  }

  //reset the chain array, should have a better way than loop over 3200 cells  
  hh_num_hits= 0; //global
  hh_num_xtrs= 0;
  good_altro_data = FALSE;  //global   
  for (j=0; j<NAL_SAMP; j++)
  {/* clear the data array */
    for (i=0; i<NUM_PADS; i++)   pad_dat[j][i]= 0.0; //global
    for (i=0; i<NUM_XTRA; i++) 	 xtra_dat[j][i]= 0.0; //global     
  }
  
  readflag = TPC_BOS_ReadTPCData();

  if( readflag <= 0 )
  {
    good_altro_data = FALSE;
    /*printf("No TPC Bank Info in this Event\n");*/
    num_bad_altro_data++;
  }
  else
  {
    //TPC data found, now start to fill up Nate's array
    //pad_dat[NAL_SAMP][NUM_PADS] and xtra_dat[NAL_SAMP][NUM_XTRA] here
    good_altro_data = TRUE;

    //read TGBI for trigger bit
    TGBI = (clasTGBI_t*) getBank(&bcs_, "TGBI");
    if(TGBI) clas_trig_type=(int)(TGBI->tgbi[0].latch1);
  }

  return readflag;
}

int TPC_BOS_ReadTPCData()
{
  int ii,jj,chan,vtdc,bankNum,ntpc,nsamples;
  int tic,pad_num,altro_signal;
  int ntoomanyhits=0;
  int imRow;
  int imCol;
  clasTPC6_t *BOSdata;
  clasTPC_t *TPC;
 
  ntpc = 0;
  hh_num_hits = 0;
 

  for(bankNum=0; bankNum<4; bankNum++)
  {
    if((BOSdata = (clasTPC6_t *) getGroup(&bcs_,"TPC6",bankNum)) > 0 )
    {
      imCol = BOSdata->bank.ncol-1;
      imRow = BOSdata->bank.nrow-1;
 
      ii=0;
      while(ii<imRow)
      {
        if (imRow < ii+2)
        {
          fprintf(stderr,"gem/read_BOS.c:  Truncated TPC6 Bank. Ignoring Event #%d\n",clas_global_eid);
          return 0;
        }

        chan = BOSdata->tpc6[ii++].data;
        vtdc = BOSdata->tpc6[ii++].data;
        nsamples = BOSdata->tpc6[ii++].data;

        if (imRow < ii+nsamples-1)
        {
          fprintf(stderr,"gem/read_BOS.c:  Truncated TPC6 Bank. Ignoring Event #%d\n",clas_global_eid);
          return 0;
        }
        
        for(jj=0; jj<nsamples; jj++)
        {
          pad_num = chan;
          tic = vtdc-jj;
          altro_signal = BOSdata->tpc6[ii++].data;

          pad_dat[tic][pad_num] = (float)altro_signal;
          if(hh_num_hits < HH_MAX_NUM_HITS)
          {
            hh_hitlist[hh_num_hits].pad = pad_num;
            hh_hitlist[hh_num_hits].t   = tic;
            hh_hitlist[hh_num_hits].q   = (float)altro_signal;
            hh_hitlist[hh_num_hits].qraw= altro_signal;
            hh_hitlist[hh_num_hits].status=HUNTCHD;
            if( (float)altro_signal < THR_Q_LINK)
              hh_hitlist[hh_num_hits].status |= HSMALLQ; 
            hh_num_hits++;
          }
          else
          {
            errors[6]++;
            ntoomanyhits++;
          }	
        }
      }
    }
    else if ((TPC = (clasTPC_t*) getGroup(&bcs_,"TPC ",bankNum)) > 0 )
    {
      imCol = TPC->bank.ncol-1;
      imRow = TPC->bank.nrow-1;
   
      for(ii=0; ii<imRow; ii++)
      {
        // note, TPC6 (above) has a check for truncation, TPC does not
        pad_num = (int)TPC->tpc[ii].id;
        tic = (int)TPC->tpc[ii].tdc/100;
        altro_signal = (int)TPC->tpc[ii].adc;
        pad_dat[tic][pad_num] = (float)altro_signal;
        if(hh_num_hits < HH_MAX_NUM_HITS)
        {
          hh_hitlist[hh_num_hits].pad = pad_num;
          hh_hitlist[hh_num_hits].t   = tic;
          hh_hitlist[hh_num_hits].q   = (float)altro_signal;
          hh_hitlist[hh_num_hits].qraw= altro_signal;
          hh_hitlist[hh_num_hits].status=HUNTCHD;
          if( (float)altro_signal < THR_Q_LINK)
            hh_hitlist[hh_num_hits].status |= HSMALLQ; 
          hh_num_hits++;
        }
        else 
        {
            errors[6]++;
            ntoomanyhits++;
        }
      }
    }
  }

  if (ntoomanyhits > 0) 
  {
    nevt_toomanyhits++;
//    printf("gem/read_BOS:  truncating data, "
//      "discarding %d hits in event #%d\n",ntoomanyhits,clas_global_eid);
  }

  return hh_num_hits;
}

//process one tpc bank ,check the hits then fill the global array
int fill_hits(clasTPC_t* pTPC, int bankindex) 
{   
  int jj;
  int altro_channel=0,altro_time=0,altro_signal=0,tic=0,pad_num=0,xtr_num=0;

  if(!pTPC) return -1; //this bank is empty, return

  // loop over this TPC bank
  for(jj=pTPC->bank.nrow-1;jj>=0;jj--)
  {
    pad_num= -999;
    //fill from the end to the start will make the data in the ascending order
    //if real data, order by chan id
    //if simulation data, order by tdc
    altro_channel = (int)pTPC->tpc[jj].id;
    altro_time    = (int)pTPC->tpc[jj].tdc;
    altro_signal  = (int)pTPC->tpc[jj].adc;
#if defined (READBOS_DEBUG) &&  (READBOS_DEBUG>=2)
      printf("TPC%d[%3d]: ID=%4d  TDC=%5d  ADC=%4d\n",
	     bankindex,jj,altro_channel,altro_time,altro_signal);
#endif
    tic =(int) altro_time/DAQ_CONVERT;
    if(altro_channel >= NAL_CHAN || altro_channel < 0) 
    {
#if defined (READBOS_DEBUG)
      printf("Row %4d: channel id from BOS_BANK is out of range: %d",jj,altro_channel);
#endif
      good_altro_data = FALSE;
      errors[7]++;
    }
    else if(altro_time%DAQ_CONVERT)
    {
#if defined (READBOS_DEBUG) &&  (READBOS_DEBUG>=1)
      printf("Row %4d: ALTRO time value is not allowed: %d",jj,altro_time);    
#endif
      good_altro_data = FALSE;
      errors[8]++;
    }
    else if( (tic >= NAL_SAMP) || (tic < 0) )
    {
#if defined (READBOS_DEBUG) &&  (READBOS_DEBUG>=1)
      printf("Row %4d: ALTRO time value is out of range: %d",jj,altro_time);
#endif
      good_altro_data = FALSE;
      errors[9]++;
    }
    else
    {
      pad_num= altro_channel; 

      if(pad_num >=0) 
      {
	pad_dat[tic][pad_num] = (float)altro_signal;
	if(hh_num_hits < HH_MAX_NUM_HITS)
	{
	  hh_hitlist[hh_num_hits].pad = pad_num;
	  hh_hitlist[hh_num_hits].t   = tic;
	  hh_hitlist[hh_num_hits].q   = (float)altro_signal;
	  hh_hitlist[hh_num_hits].qraw= altro_signal;
	  hh_hitlist[hh_num_hits].status=HUNTCHD;
	  if( (float)altro_signal < THR_Q_LINK) 
	    hh_hitlist[hh_num_hits].status |= HSMALLQ; /* charge too small */
	  hh_num_hits++;
	}
	else errors[6]++;
      }
      else if(xtr_num >=0)
      {
	xtra_dat[tic][xtr_num]= (float)altro_signal;
	if(hh_num_xtrs < HH_MAX_NUM_XTRS)
	{
	  hh_xtrlist[hh_num_xtrs].pad= xtr_num;
	  hh_xtrlist[hh_num_xtrs].t= tic;
	  hh_xtrlist[hh_num_xtrs].q= (float)altro_signal;
	  hh_xtrlist[hh_num_hits].qraw= altro_signal;
	  hh_num_xtrs++;
	}
	else errors[10]++;
      }
    }//end of filling one hit
  } //end of TPC bank
  return 1;
}
