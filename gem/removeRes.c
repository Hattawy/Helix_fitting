#include <stdio.h>
#include "rtpc.h"

int removeRes(int tracknum, int chainnum)
{
  /////////////////////////////////////////////////////////////////////////
  //this routine calculates the residual of each point, removes the point
  //if it is too far away from the fit, then returns 1 if the average of the 
  //residuals is > RES_CUT and 0 if it isn't.  (WARNING - actually returns
  //#hits on the chain after rejecting ones with large residuals).
  ////////////////////////////////////////////////////////////////////////


  #define KILL_NUM 4
  
  typedef struct {
    int hitnum;
    float del;
  } out_list;
  
  out_list kill[KILL_NUM];
  int hold_hitnum;
  float hold_del;
  float residual, resSum, resAvg;
  int jj, ii, kk, n;
  int decideFate;
  int removeLink;
  int num_links_to_remove;
  int KILLDEBUG = FALSE;

  n = 0;
  num_links_to_remove = 0;
  resSum = 0.0;
  removeLink = FALSE;
   
  if(KILLDEBUG) printf("starting removeRes with chainnum = %d, hh_num_hits = %d, num_hits_this_chain = %d\n",
		    chainnum, hh_num_hits, num_hits_this_chain[chainnum]);

  if(KILLDEBUG)
  {
    for(ii = 0; ii<num_hits_this_chain[chainnum]; ii++)
     printf("chain %d point %d = %d\n", chainnum, ii, chain_hits[chainnum][ii]);
  }
  //initialize the kill list
  for (ii=0; ii<KILL_NUM; ii++)
  {
    kill[ii].hitnum = -1000;
    kill[ii].del = -1000;
  }
  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {
    if(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)
    {
      residual = hh_hitlist[chain_hits[chainnum][jj]].del;
      n++;
      if(KILLDEBUG) printf("hit num %d on chain %d has residual %f\n",jj, chainnum, residual);
      resSum += residual;
      //if this link is > RESCUT, see if it is one of the top four largest, if it is, add
      //it to the kill list
      if(residual > RESCUT) 
      {
        decideFate = TRUE;
        ii=0;
        while(decideFate)
    	{
    	  //exit loop if we are at the end of the kill list
    	  if(ii==KILL_NUM) decideFate = FALSE;
    	  //now check to see if the point has a residual larger than 
    	  //the current point on the kill list
    	  else if(residual > kill[ii].del)
    	  {
    	    hold_hitnum = kill[ii].hitnum;
    	    removeLink = TRUE;
    	    hold_del = kill[ii].del;
    	    kill[ii].del = residual;
    	    kill[ii].hitnum = chain_hits[chainnum][jj];
    	    decideFate = FALSE; //we want to exit the loop if we found a value
                                //which is less than the current del

    	    //now drop the last point from the kill list, then traverse
    	    //the list backwards till we come to the point just after
    	    //the point we added
    	    for(kk=KILL_NUM-1; kk>ii; kk--)
    	    {
    	      if(kk == ii+1)
    	      {
    	        kill[kk].hitnum = hold_hitnum;
    	        kill[kk].del    = hold_del;
    	      }
    	      else
    	      {
    	        kill[kk].hitnum = kill[kk-1].hitnum;
    	        kill[kk].del    = kill[kk-1].del;
    	      }
    	    }
    	  }
    	  else ii++;
    	}//end of while loop over kill list
      }//end of check to see if this is an outlier
    }//end of status bit check
  }//end of loop over hits in the chain
  //what we should have now is a list of points held in the array kill which
  //is ordered [0, KILL_NUM-1] from largest residual to smallest over RESCUT
  
  if(removeLink)
  {
    for(ii=0; ii<KILL_NUM; ii++)
    {
      if(kill[ii].hitnum != -1000)
      {
        hh_hitlist[kill[ii].hitnum].status |= HREMOVD;
        hh_hitlist[kill[ii].hitnum].status &= ~HISUSED;
        num_links_to_remove++;
        if(KILLDEBUG) printf("removed link %d from chain %d\n",
      		  kill[ii].hitnum, chainnum);
      }
    }
  }

  if(KILLDEBUG)
    for(ii = 0; ii<num_hits_this_chain[chainnum]; ii++)
      if(hh_hitlist[chain_hits[chainnum][ii]].status & HISUSED)
        printf("chain %d point %d = %d\n", chainnum, ii, chain_hits[chainnum][ii]);
  if(KILLDEBUG) printf("after removing %d outliers, chain %d has npts = %d\n",
		      num_links_to_remove,chainnum,num_hits_this_chain[chainnum] - num_links_to_remove);
  resAvg = resSum/(float)n;


  return (num_hits_this_chain[chainnum] - num_links_to_remove);
}
