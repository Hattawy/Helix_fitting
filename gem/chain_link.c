#include <stdio.h>
#include <math.h>
#include "rtpc.h"

// nab modified 2012
// - no more checking data size, done now in gem_evnt
// - no more fortran calls (they were rounding floats to ints)
// - no more truncating when one chain has more than MAX_HITS_ON_CHAIN
//   - now just ignore that chain (and all its hits for future chains)

int chain_link()
{
  /* find tracks by stringing together links of hits that are close */

  /*
     Hits are stored in a structure that can be used to point into rawCYL:		
     hh_num_hits = number of hits in the RTPC for this event
     chain_hits[chain no.][hit# on chain] contains pointers to the hits in hh_hitlist on each chain
     num_chains = number of chaings found
     num_hits_this_chain[chain no.] = number of hits on each chain
     int   hh_hitlist[hh_num_hits].pad= pad_num -- the pad number (0-3199)
     int   hh_hitlist[hh_num_hits].t  = tic     -- the time of the hit (0-99)
     int   hh_hitlist[hh_num_hits].status       -- the unused/used status of this hit. Initially 0.
     HITUNAV means it is not available for linking.
     float hh_hitlist[hh_num_hits].z            -- z position of the pad center (mm)
     float hh_hitlist[hh_num_hits].r            -- radial position of the hit (mm)
     float hh_hitlist[hh_num_hits].phi          -- phi of the hit (after Lorentz correction)
     float hh_hitlist[hh_num_hits].q            -- signal pulse height for this hit
     */
  int ii,jj,anchor_hit, seed_hit, next_hit, seed_index;
  float pseed[3], pnext[3]; /* place to put position 3-vectors of points */
  float separation;
  int chain_debug=FALSE;//TRUE;//FALSE;

  num_chains=0;

  /* Print out what we started with */
  if(chain_debug)
  {
    printf("***************************** \n");
    printf(" Entering chain_link with %d hits (eid = %d)\n",hh_num_hits,clas_global_eid);
    printf(" HitNum   Pad  Time Sta     r      z     phi      q\n");
    for(jj=0; jj<hh_num_hits; jj++)
    {
      if (jj >= HH_MAX_NUM_HITS) break;
      printf(" %6d %5d %5d %3d %6.1f %6.1f %6.3f %6.0f \n",              
          jj,
          hh_hitlist[jj].pad,
          hh_hitlist[jj].t,
          hh_hitlist[jj].status,
          hh_hitlist[jj].r,
          hh_hitlist[jj].z,
          hh_hitlist[jj].phi,
          hh_hitlist[jj].q   );
    }
  }


  for(anchor_hit=0; anchor_hit < hh_num_hits; anchor_hit++)
  {

    if ((hh_hitlist[anchor_hit].status & HITUNAV))
    {
      if(chain_debug) 
      {
        printf("Hit %d unavailable as anchor: status=%d mask=%d\n",
            anchor_hit,hh_hitlist[anchor_hit].status,HITUNAV);
      }
      continue;
    }
    if( num_chains >= MAX_NUM_CHAINS)
    {
      fprintf(stderr,"gem/chain_link: too many chains, truncating event#%d\n",clas_global_eid);
      errors[4]++;
      return -1;
    }

    num_hits_this_chain[num_chains]=1;
    chain_hits[num_chains][0]= anchor_hit;

    if(chain_debug) printf("Trying new chain (number %d) anchored at hit %d\n",num_chains,anchor_hit);

    hh_hitlist[anchor_hit].status |= HISUSED; /* at least temporarily, mark the anchor hit as used */

    for(seed_hit=0; seed_hit < num_hits_this_chain[num_chains];seed_hit++)
    {

      if(chain_debug) printf(" using seed hit %d\n",chain_hits[num_chains][seed_hit]);

      /* get (x,y,z) of the seed hit*/
      seed_index= chain_hits[num_chains][seed_hit];
      pseed[0]= hh_hitlist[seed_index].r * cos(hh_hitlist[seed_index].phi);
      pseed[1]= hh_hitlist[seed_index].r * sin(hh_hitlist[seed_index].phi);
      pseed[2]= hh_hitlist[seed_index].z;

      /* check this seed hit against ALL unused hits */
      for(next_hit= 0 ; next_hit < hh_num_hits; next_hit++)
      {

        if ( hh_hitlist[next_hit].status & HITUNAV) continue;

        /* If hit not on the same RTPC side */
        if ( (hh_hitlist[next_hit].pad   -1600)*
            (hh_hitlist[seed_index].pad -1600) < 0 ) continue;

        pnext[0]= hh_hitlist[next_hit].r * cos(hh_hitlist[next_hit].phi);
        pnext[1]= hh_hitlist[next_hit].r * sin(hh_hitlist[next_hit].phi);
        pnext[2]= hh_hitlist[next_hit].z;

        // used to be done with nvmod and nvsub fortran routines:
        separation = sqrt( pow(pseed[0]-pnext[0],2) + 
                           pow(pseed[1]-pnext[1],2) +
                           pow(pseed[2]-pnext[2],2) );

        //      if(chain_debug) printf("... separation is %5.1f\n",separation);

        if(separation >= MAX_LINK_SEP) continue;

        // next_hit is close enough to seed_hit, so try to add it to the chain:
        
        if(chain_debug) printf("...adding hit %d to this chain\n",next_hit);

//        if (num_hits_this_chain[num_chains] >= MAX_HITS_ON_CHAIN)
//        {
//          return -1;
//        }

        // mark it as used:
        hh_hitlist[next_hit].status |= HISUSED;
       
        // If there's too many hits on the chain, just leave them marked as used
        // so future chains can't reuse it, but don't add it to the chain array.
        // This entire chain will be ignored later.
        if (num_hits_this_chain[num_chains] >= MAX_HITS_ON_CHAIN) continue;

        // add hit to array:
        chain_hits[num_chains][num_hits_this_chain[num_chains]]=next_hit;
        num_hits_this_chain[num_chains]++;

      }
    }
     
    if( num_hits_this_chain[num_chains] > 1) 
    {
   
      // If there's too many hits on the chain, ignore it:
      if (num_hits_this_chain[num_chains] >= MAX_HITS_ON_CHAIN) continue;
      
      if(chain_debug)
      {
        printf("....KEEPING THIS CHAIN #%d. %d hits \n",
            num_chains,num_hits_this_chain[num_chains]);
        for (jj=0; jj<num_hits_this_chain[num_chains]; jj++)
          printf(" %d", chain_hits[num_chains][jj]);
        printf("\n\n");
      }
      
      num_chains++; /* save chain if >1 hit on it */
    
    }
    else
    {
        hh_hitlist[anchor_hit].status &= ~HISUSED;
    } /* no chain -> release the anchor hit */
  }
  if(GEMDEBUG)
  {
    /* print out a summary of found chains */
    for (ii=0; ii<num_chains; ii++)
    {
      printf("\n");
      printf("Hits used on chain %d (%d):",ii,num_hits_this_chain[ii]);
      for (jj=0; jj<num_hits_this_chain[ii]; jj++)
        printf(" %d", chain_hits[ii][jj]);
    }
    printf("\n");
  }
  return num_chains;
}
