#include <stdio.h>
#include <math.h>
#include "rtpc.h"

// nab modified 2012
// - returning negative for bad tracks immediately.
// - not finding fiducial entry/exit points, weren't used anyway.
// - stepping now starts at smin and ends at smax since were're
//   not finding fiducial entry/exit anymore, and vtl (path length
//   not correceted for badpads) is and was just smax-smin.
// - step size decreased, but no slower due to above changes.
// - implemented live_length (path_length corrected for bad pads).
//
// nab modified 2014
// - add truncated path length calculation (live_lengthL)
//

float track_charge(int tracknum,float *vtl,float *sdist,float *edist)
/* for track number 'tracknum' find the first and last hits along the helix 
 * already fitted to it.  Returned  value  is the sum of hits'  pulse-heights  
 * along the helix within the  visible  portion of the rtpc. Also return:
	*vtl= length along helix from first hit to last  hit
	*sdist= distance from rtpc entrance point  to the first  hit
	*edist= distance, along helix, from  rtpc exit point to last hit

	'first' and  'last'  hits are defined  in the sense  that  the helix 
        originates at the fitted point close[tracknum].  */
{
  float dummy, sumq, sig, altro_time_at_pad, prow, pcol;
  int num_x, pad_num, jj;
  cylindrical point, atpad;

  // step size 0.25 was the original (nab)
  const float sig_step=0.1;

  int nsteps=0;
  trackv[tracknum].live_length=0.0;
  trackv[tracknum].live_lengthL=0.0;
  trackv[tracknum].live_lengthA=0.0;
  trackv[tracknum].live_lengthAL=0.0;
  trackv[tracknum].bad=0.0;

  // check if point of closest approach is inside tracking volume:
  if(inside_rtpc(x_close[tracknum],y_close[tracknum],z_close[tracknum]))
  {
    // This never happens for EG6 (we have beamline helix constraint)
    if(!findwall(tracknum,-1, &trackv[tracknum].sig_entr,&dummy)) errors[14]++;
    return -1.;
  }
  
  // get boundary points of this track: where it enters/leaves rtpc:
  num_x=findwall(tracknum, 1, &trackv[tracknum].sig_entr, &trackv[tracknum].sig_exit);

  if(num_x<1) 
  {
    // This happens for about 0.2% of EG6 "tracks".
    errors[13]++;
    return -2.;
  }
  // num_x=1 for track with small radius of curvature whose
  // helix intersects the cathode but no other boundary.

  // calculate track charge (2014, updated to also do truncated versions):
  sumq = sumtrackq(tracknum, &trackv[tracknum].smin,   &trackv[tracknum].smax,
  		             &trackv[tracknum].t_smin, &trackv[tracknum].t_smax);


  // points of first and last detected charge on the track:
  sig2xyz(tracknum, trackv[tracknum].smin, &point_appear);
  sig2xyz(tracknum, trackv[tracknum].smax, &point_disappear);
  point_appear.t    = trackv[tracknum].t_smin;
  point_disappear.t = trackv[tracknum].t_smax;
  
  // visible track length:
  *vtl  = trackv[tracknum].smax     - trackv[tracknum].smin;
  // missing track length: from entry point to first signal:
  *sdist= trackv[tracknum].smin     - trackv[tracknum].sig_entr;
  // disappearing track length: from last signal to exit point:
  *edist= trackv[tracknum].sig_exit - trackv[tracknum].smax;

  // stop wasting time:
  if (fabs(*sdist)>50 || fabs(*edist)>50) return -3.;


  // load arrays of indices of hits used in dedx truncations:
  const int chainnum=track_chain[tracknum];
  int nindA=0,nindAL=0,nind=0;
  int indA[MAX_HITS_ON_CHAIN],indAL[MAX_HITS_ON_CHAIN];

//  int nindL=0;
//  int indL[MAX_HITS_ON_CHAIN];

  for (jj=0; jj< num_hits_this_chain[chainnum]; jj++)
  {
    if( !(hh_hitlist[chain_hits[chainnum][jj]].status & HISUSED)) continue;

    nind++;

    if (hh_hitlist[chain_hits[chainnum][jj]].status & HTRUNCA)
      indA[nindA++] = chain_hits[chainnum][jj];

    if (hh_hitlist[chain_hits[chainnum][jj]].status & HTRUNCAL)
      indAL[nindAL++] = chain_hits[chainnum][jj];

//    if (hh_hitlist[chain_hits[chainnum][jj]].status & HTRUNCL)
//      indL[nindL++] = chain_hits[chainnum][jj];
  }


/*
     fprintf(stderr,"%5d %5d %5d %5d %12.4f %12.4f %12.4f\n",
     nind,nindL,nindA,nindAL,
          (float)nindL/nind,
          (float)nindA/nind,(float)nindAL/nind);
*/

/*
  fprintf(stderr,"%d %d %d %d %d\n",tracknum,nind,nindL,nindA,nindAL);
  for (jj=0; jj<nindL; jj++)
    fprintf(stderr,"%d ",indL[jj]);
  fprintf(stderr,"\n");
  for (jj=0; jj<nindA; jj++)
    fprintf(stderr,"%d ",indA[jj]);
  fprintf(stderr,"\n");
  for (jj=0; jj<nindAL; jj++)
    fprintf(stderr,"%d ",indAL[jj]);
  fprintf(stderr,"\n");
*/

  // calculate active path length:
  for (sig=trackv[tracknum].smin; sig<=trackv[tracknum].smax; sig+=sig_step)
  {
    // get point on helix
    sig2rpz(tracknum, sig, &point);

    if (!inside_rtpc_cyl(point)) continue;

    // get corresponding point on gem surface
    if (proj_to_pads_fit_eg6(point, &atpad, &altro_time_at_pad)) continue;
    
    pad_num= prowcol(atpad.z, atpad.phi, &prow, &pcol);
  
    if (pad_num>=0)
    {
      nsteps++;
      if (dead_pad(pad_num)) 
      {
          trackv[tracknum].bad += 1.;
      }
      else
      {
        // Non-truncated path-length:
        trackv[tracknum].live_length += sig_step;

        if (sig >= trackv[tracknum].sminL &&
            sig <= trackv[tracknum].smaxL)
        {
          // Path-length-truncated path length:
          trackv[tracknum].live_lengthL += sig_step;

          // Path-length- and ADC- truncated path length
          // (test for hit matching, should binary search for speed)
          for (jj=0; jj<nindAL; jj++)
          {
            if (hh_hitlist[indAL[jj]].pad == pad_num &&
                hh_hitlist[indAL[jj]].t   == altro_time_at_pad)
            {
              trackv[tracknum].live_lengthAL += sig_step;
              break;
            }
          }


        }

        // ADC-truncated path length
        // (test for hit matching, should binary search for speed)
        for (jj=0; jj<nindA; jj++)
        {
          if (hh_hitlist[indA[jj]].pad == pad_num &&
              hh_hitlist[indA[jj]].t   == altro_time_at_pad)
          {
            trackv[tracknum].live_lengthA += sig_step;
            break;
          }
        }


      }
    }
  
  }



  // divide # badpad integration steps by total number of steps:
  trackv[tracknum].bad = nsteps>0 ? trackv[tracknum].bad/nsteps : 1.;

  return(sumq);




  // getting rid of the rest of this.
  // everything that was actually used anywhere else in the code is above here.
  // (nab 2012)





  /* The physical RTPC volume defined by findwall.c is not the same as the 
     sensitive volume of the chamber because of the Lorentz angle. Here we
     are going to move along the helix between the physical entry and exit
     points and find the first and last points on the helix that are within
     the actual sensitive volume. The method is to find the pad-number to
     which electrons would drift. A valid pad number implies that the point
     is indeed within the sensitive volume, since the pad-finding code knows
     about the electron drift-paths. 

     While we are stepping along the helix, we are going to add up the 
     accumulated step distances that would project electrons onto pads that have
     NOT been identified as dead. This gives us a measure of the "live" track
     length.
  */

/*
  //float start_sig, stop_sig;
  //float xxp, yyp, zzp;
  //cartesian point_fid_entr, point_fid_exit;
  //int found_fid_entr=FALSE;
  //int found_fid_exit=FALSE;
  //int nFalseExit=0;
  //int nFalseEntrance=0;
  //int inside=0;
  
  start_sig = trackv[tracknum].sig_entr - 10.0;
  stop_sig  = trackv[tracknum].sig_exit + 10.0;

  trackv[tracknum].sig_fid_entr = trackv[tracknum].sig_entr;
  trackv[tracknum].sig_fid_exit = trackv[tracknum].sig_exit;  

  if (stop_sig <= start_sig)
  {
    //fprintf(stderr,"stop_sig<=start_sig %d %d %f\n",clas_global_eid,tracknum,r_0[tracknum]);
    return -1.;
  }

  for (sig = start_sig ; sig < (stop_sig+sig_step) ;  sig += sig_step)
  {
    if (found_fid_exit) break;
    if (sig-start_sig >= 500.0) break;

    // get (r,phi,z) point along helix:
    sig2rpz(tracknum, sig, &point);

    // get (r,phi,z) on gem surface: (reverse drift path calibration)
    proj_to_pads_error = proj_to_pads_fit_eg6(point, &atpad, &altro_time_at_pad);

    if(proj_to_pads_error)
    {
      if(GEMDEBUG) printf("proj_to_pads_error, killing track");
      // we don't want to break here because we are starting before drift region
      // where proj_to_pads' calibrations will fail, just wait till we reach valid radius
      // continue;
      break;
    }

    // (x,y,z) point along helix:
    xxp= point.s*cos(point.phi);
    yyp= point.s*sin(point.phi);
    zzp= point.z;

    inside = inside_rtpc(xxp,yyp,zzp);

    // get pad# from (z,phi) point on gem surface:
    pad_num= prowcol(atpad.z, atpad.phi, &prow, &pcol);

    //  sum up track length in actual live volume of rtpc
    if (pad_num>=0 && inside)
    {
      // like the other path length calculation (vtl), ignore the
      // missing portions at beginning and end of track
      if (sig>=trackv[tracknum].smin && sig<=trackv[tracknum].smax)
      {
        nsteps++;
        if (dead_pad(pad_num)) trackv[tracknum].bad += 1.;
        else                   trackv[tracknum].live_length += sig_step;
      }
    }




    if(!found_fid_entr)
    {
      if(GEMDEBUG) printf("looking for entrance %d: sig= %6.1f, (s,phi,z)= (%6.1f,%6.3f,%6.1f) ",nFalseEntrance,sig, point.s, point.phi, point.z);
      if( pad_num>=0 && inside )
      {
        found_fid_entr= TRUE;
        trackv[tracknum].sig_fid_entr= sig;
      }
      else nFalseEntrance++;
      if(nFalseEntrance > 100000)
      {
        if(GEMDEBUG) printf("Tried more than 100000 to find an entrance, quitting");
        break;
      }
    }
    else
    {
      if(GEMDEBUG) printf("looking for exit: sig= %6.1f, (s,phi,z)= (%6.1f,%6.3f,%6.1f) ",sig, point.s, point.phi, point.z); 
      if( !( pad_num>=0 && inside ) ) // Exited sensitive volume
      {
        found_fid_exit= TRUE;
        trackv[tracknum].sig_fid_exit= sig-sig_step;
      }
      if(GEMDEBUG) printf(" .... found_exit= %d\n",found_fid_exit); 
      if(!found_fid_exit) nFalseExit++;
      if(nFalseExit > 100000) 
      {
        if(GEMDEBUG) printf("Tried more that 100000 to find an exit, quitting");
        break;
      }
    }
  }


  if (!found_fid_exit)
  { 
    if(GEMDEBUG)
    {
      printf(">>>WARNING: No fid-exit point for track# %d in event# %d ",tracknum,clas_global_eid); 
      printf("between sig= %6.1f and %6.1f -- %f\n",start_sig, stop_sig, r_0[tracknum]); 
    }
    helix_error[tracknum] = 1;
  }
  if (!found_fid_entr)
  { 
    if(GEMDEBUG)
    {
      printf(">>>WARNING: No fid-entry point for track# %d in event# %d ",tracknum,clas_global_eid); 
      printf("between sig= %6.1f and %6.1f -- %f\n",start_sig, stop_sig,r_0[tracknum]); 
    }
    helix_error[tracknum] = 2;
  } 


  // remember noteworthy points along this track
  sig2xyz(tracknum, trackv[tracknum].sig_entr,     &point_entry);
  sig2xyz(tracknum, trackv[tracknum].sig_exit,     &point_exit);
  sig2xyz(tracknum, trackv[tracknum].sig_fid_entr, &point_fid_entr);
  sig2xyz(tracknum, trackv[tracknum].sig_fid_exit, &point_fid_exit);
  sig2xyz(tracknum, trackv[tracknum].smin,         &point_appear);
  sig2xyz(tracknum, trackv[tracknum].smax,         &point_disappear);
  point_appear.t    = trackv[tracknum].t_smin;
  point_disappear.t = trackv[tracknum].t_smax;

  *vtl  = trackv[tracknum].smax     - trackv[tracknum].smin; // visible track length
  *sdist= trackv[tracknum].smin     - trackv[tracknum].sig_entr; // missing track length: from entry point to first signal
  *edist= trackv[tracknum].sig_exit - trackv[tracknum].smax; // disappearing track length: from last signal to exit point




  if(GEMDEBUG)
  {
    printf("Track enters rtpc volume at s=%5.1f (x,y,z)=(%6.1f,%6.1f,%6.1f)\n",
        trackv[tracknum].sig_entr,point_entry.x,point_entry.y,point_entry.z);

    printf(" enters sensitive volume at s=%5.1f (x,y,z)=(%6.1f,%6.1f,%6.1f)\n",
        trackv[tracknum].sig_fid_entr,    point_fid_entr.x,point_fid_entr.y,point_fid_entr.z);

    printf("         becomes visible at s=%5.1f (x,y,z)=(%6.1f,%6.1f,%6.1f)\n",
        trackv[tracknum].smin,    point_appear.x,point_appear.y,point_appear.z);

    printf("              disappears at s=%5.1f (x,y,z)=(%6.1f,%6.1f,%6.1f)\n",
        trackv[tracknum].smax,    point_disappear.x,point_disappear.y,point_disappear.z);

    printf(" leaves sensitive volume at s=%5.1f (x,y,z)=(%6.1f,%6.1f,%6.1f)\n",
        trackv[tracknum].sig_fid_exit,    point_fid_exit.x,point_fid_exit.y,point_fid_exit.z);

    printf("  and leaves rtpc volume at s=%5.1f (x,y,z)=(%6.1f,%6.1f,%6.1f)\n",
        trackv[tracknum].sig_exit,point_exit.x,point_exit.y,point_exit.z);

    printf("Track runs from s=%5.2f to %5.2f (%5.1fmm)\n",
        trackv[tracknum].sig_entr,trackv[tracknum].sig_exit,
        trackv[tracknum].sig_exit-trackv[tracknum].sig_entr);
    printf("      missing length= %5.2f disappearing length= %5.2f charge=%6.1f\n",*sdist,*edist,sumq);
  }

  return(sumq);
*/

}
