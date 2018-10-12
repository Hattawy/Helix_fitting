#include <math.h>
#include <stdio.h>
#include <ntypes.h>
#include <bostypes.h>
#include <bosddl.h>
#include "rtpc.h"

// nab modified 2012, 2014

extern int clas_trig_type;

int gem_evnt_()
{
  int nRTPC,iGCPB,RTPC_pids[5];
  float RTPC_moms[5];
  gcpb_t gcpb[EVT_TRKS];
  rtpc_t rtpc[EVT_TRKS];

  int ii,jj,kk;
  int npads_hit, npads_track;
  int original_num_pts,npts_tofit;
  float vtl;            // visual track length
  float sdist;          // distance from start of track to cathode
  float edist;          // distance from end of track to some chamber boundary
  float b_z;            // magnetic field value from map
  float visible_charge; // total charge for a track
  float p_perp, p_para, p_x, p_y, p_tot, sig0, qtot, dx;
  cartesian vertex; // point on helix closest to beamline (not necessarily same as z axis

  tpc_hits_.tpc_hits_nb=0; // number of rows in TPCC ntuple
  
  // put data into array hh_hitlist.pad,t,q,status:
  if (read_BOS() <= 0) return 0;

  if (hh_num_hits > HH_MAX_NUM_HITS)
  {   // impossible
      fprintf(stderr,"gem/gem_evnt:  BOS INPUT ERROR.  Skipping Event.\n");
      return 0;

  }

  num_events++;
  
  gem_print_stats();
  
  hits_to_space(); // drift paths calibration (set hh_hitlist.z,r,phi)

  scale_hits();    // gain calibration (scale hh_hitlist.q)
  
  clean_hits();    // set hh_hitlist.status (as of 2014, must be called AFTER scale_hits)
  
  chain_link();    // fill chain_hits with pointers to hh_hitlist
  
  sum_hit_charges(&npads_hit, &qtot);  // sum adcs for the whole event

  ntracks=0;
  for(ii=0; ii<num_chains; ii++) 
  {
    if (num_hits_this_chain[ii] >= MIN_LINKS_ON_TRACK) ntracks++;
  }
  if (ntracks == 0) return 0;


  //  anachains();


  ntracks=0;  // number of rows in GCPB bank


  for(ii=0; ii<num_chains; ii++)
  {

    if(ntracks >= EVT_TRKS)
    {
      ntracks = EVT_TRKS;
      printf("gem/gem_evnt:  Truncating to %d tracks in event# %d.\n", EVT_TRKS,clas_global_eid);
      break;
    }

    if (num_hits_this_chain[ii] < MIN_LINKS_ON_TRACK) continue;

    // store index pointer from this track# to chain_hits array index:
    track_chain[ntracks]= ii;

    // run helix fit:
    original_num_pts = num_hits_this_chain[track_chain[ntracks]];
    helix_fit(fit_track_to_beamline, ntracks, track_chain[ntracks]);

    // remove large residuals, rerun fit:
    npts_tofit = removeRes(ntracks,track_chain[ntracks]);
    if (npts_tofit < MIN_LINKS_ON_TRACK) continue;
    if (original_num_pts != npts_tofit)
    {
      if (GEMDEBUG) printf("we removed outliers so...second call to helix fit\n");
      helix_fit(fit_track_to_beamline, ntracks, track_chain[ntracks]);
    }
    
    // bad reult from helix_fit:
    if (isnan(chi2[ntracks]) != 0) continue;
    if (chi2[ntracks]>30. || fabs(r_0[ntracks]) < 1.) continue;
    
    npads_track = countPads(ntracks,track_chain[ntracks]);
    if (npads_track < 4) continue;

    // charge/path_lengths analysis:
    visible_charge = track_charge(ntracks, &vtl, &sdist, &edist);
    
    // bad result from track_charge:
    if (visible_charge < 1.e-3) continue;
    if (isnan(vtl) != 0 || isnan(edist) != 0 || isnan(sdist) != 0) continue;
    if (fabs(edist) > 50 || fabs(sdist) > 50) continue;


    track_phi= phi0[ntracks];
    track_theta = PIBY2 - atan(dzds[ntracks]); /* atan returns -PI/2 to PI/2 */

    b_z = bofz((point_appear.z + point_disappear.z)/2.0); //gives b_z in kiloguass
    b_z /= 10.0; //converts to Tesla
    
    if (b_z < 0.0) p_perp = (0.3)*(3.876)*fabs(r_0[ntracks]); 
    else           p_perp = (0.3)*(float)b_z*fabs(r_0[ntracks]); 
    p_perp= fabs(p_perp);
    p_x = p_perp * cos(track_phi);
    p_y = p_perp * sin(track_phi);
    p_para = p_perp/tan(track_theta);
    p_tot = sqrt(p_perp*p_perp + p_para*p_para); 

    trackv[ntracks].p_tot  = p_tot;
    trackv[ntracks].p_corr = 0;
    trackv[ntracks].p_x    = p_x;
    trackv[ntracks].p_y    = p_y;
    trackv[ntracks].p_z    = p_para;
    trackv[ntracks].xdca   = x_close[ntracks];
    trackv[ntracks].ydca   = y_close[ntracks];
    trackv[ntracks].zdca   = z_close[ntracks];
    trackv[ntracks].sign   = r_0[ntracks]/fabs(r_0[ntracks]);
    trackv[ntracks].r_0    = r_0[ntracks];
    trackv[ntracks].theta  = track_theta;
    trackv[ntracks].phi    = track_phi;
    trackv[ntracks].dqdx   = visible_charge/vtl;
    trackv[ntracks].npts   = npts_tofit;
    trackv[ntracks].visq   = visible_charge;
    trackv[ntracks].edist  = edist;
    trackv[ntracks].sdist  = sdist;
    trackv[ntracks].vtl    = vtl;
    trackv[ntracks].chi2   = chi2[ntracks];
    trackv[ntracks].fitqual= -99999.9;


    /* The helix-fitter and the derivations above report the track 
     * parameters at the point of closest approach to the BONUS 
     * Z-AXIS. For comparison with CLAS (UserAna output), we 
     * need to know them at the point of closest approach to 
     * the CLAS Z-axis. For doing physics (e.g. scattering angle) 
     * analysis, we will want them relative to the beamline at 
     * the point of closest approach to that beamline.  Here we 
     * compute and save helical track parameters at the closest 
     * point to the beamline, but still in BONUS coordinates. */
    sig0= lineandhelix(x_0[ntracks],y_0[ntracks],z_close[ntracks],
                       r_0[ntracks],dzds[ntracks],
                       bmx0-XBONUS,bmy0-YBONUS,0.0,bmdxdz,bmdydz);

    sig2xyz(ntracks, sig0, &vertex);

    trackv[ntracks].xv=vertex.x;
    trackv[ntracks].yv=vertex.y;
    trackv[ntracks].zv=vertex.z;
    trackv[ntracks].thv=track_theta;
    trackv[ntracks].phv=phrng(atan2(vertex.y-y_0[ntracks], 
          vertex.x-x_0[ntracks])+trackv[ntracks].sign*PIBY2 );


    if(writePadUses)  // study q on each pad
      track_qpad(ntracks, p_tot, trackv[ntracks].dqdx, sdist, vtl);


    if (fabs(sdist)<10 && fabs(edist)<10 && fabs(chi2[ntracks])<10)
    {
      // fill TPCC ntuple variables:
      for(jj=0 ; jj < num_hits_this_chain[track_chain[ntracks]] ; jj++)
      {
        if (tpc_hits_.tpc_hits_nb>=1000) // array limit in ntuple
        {
          //        fprintf(stderr,"gem/gem_evnt:  TOO MANY TPCC HITS\n");
          break;
        }
        kk = chain_hits[track_chain[ntracks]][jj];
        if(hh_hitlist[kk].status & HISUSED)
        {
          tpc_hits_.chan_nb[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].pad;
          tpc_hits_.adc[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].qraw;
          tpc_hits_.cadc[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].q;
          tpc_hits_.tdc[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].t;
          tpc_hits_.dist[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].dist;
          tpc_hits_.qres[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].res_q;
          tpc_hits_.qresr[tpc_hits_.tpc_hits_nb] = hh_hitlist[kk].res_qr;
          tpc_hits_.track_nb[tpc_hits_.tpc_hits_nb] = ntracks; // pointer to GCPB
          tpc_hits_.tpc_hits_nb ++ ;
        }
      }
    }

    // print_track(ntracks);


    // store info for GCPB bank:
    gcpb[ntracks].pid        = 0;//1
    gcpb[ntracks].vert.x     = trackv[ntracks].xv;//2
    gcpb[ntracks].vert.y     = trackv[ntracks].yv;//3
    gcpb[ntracks].vert.z     = trackv[ntracks].zv;//4
    gcpb[ntracks].dedx       = visible_charge/vtl;//5
    gcpb[ntracks].p.x        = p_x;//6
    gcpb[ntracks].p.y        = p_y;//7
    gcpb[ntracks].p.z        = p_para;//8
    gcpb[ntracks].p_tot      = p_tot;//9
    gcpb[ntracks].x2         = chi2[ntracks];//10
    gcpb[ntracks].theta      = trackv[ntracks].thv;//11
    gcpb[ntracks].q          = (float)visible_charge;//12
    gcpb[ntracks].dca        = dca[ntracks];//13
    gcpb[ntracks].index      = clas_global_eid;//(int)num_events; //14
    gcpb[ntracks].phi        = trackv[ntracks].phv; //15
    gcpb[ntracks].vtl        = vtl; //16
    gcpb[ntracks].sdist      = sdist; //17
    gcpb[ntracks].edist      = edist; //18
    gcpb[ntracks].npts       = npts_tofit;//19 
    gcpb[ntracks].r_0        = r_0[ntracks]; //20
    gcpb[ntracks].fiterr     = helix_error[ntracks]; //21
    gcpb[ntracks].tothits    = hh_num_hits; //22
    gcpb[ntracks].npd_track  = npads_track;//23
    gcpb[ntracks].npd_event  = npads_hit;//24
    gcpb[ntracks].bonus_bits = 0;//bonus_bits;//25
    gcpb[ntracks].q_tot      = qtot;//26
    gcpb[ntracks].x_start    = point_appear.x;//27
    gcpb[ntracks].y_start    = point_appear.y;//28
    gcpb[ntracks].z_start    = point_appear.z;//29
    gcpb[ntracks].x_end      = point_disappear.x;//30
    gcpb[ntracks].y_end      = point_disappear.y;//31
    gcpb[ntracks].z_end      = point_disappear.z;//32


    totTrkd++; // # tracks total (global variable)
    ntracks++; // # tracks in this one event


  }//end of loop over number of chains found



  // fill the BOS banks:

  if (ntracks>0)
  {
    clasGCPB_t* GCPB = (clasGCPB_t *) 
      makeBank (&bcs_, "GCPB", 0, sizeof (gcpb_t) / 4, ntracks);

    // number of tracks to go in RTPC bank:
    nRTPC=0;

    for (iGCPB=0; iGCPB<ntracks && iGCPB<EVT_TRKS; iGCPB++)
    {
      // fill GCPB bank row:
      GCPB->gcpb[iGCPB]=gcpb[iGCPB];

      // select best tracks to go in RTPC bank:
      if (
          //gcpb[iGCPB].edist > -2      && gcpb[iGCPB].edist < 5      &&
          //gcpb[iGCPB].dca   > -0.0002 && gcpb[iGCPB].dca   < 0.0004 &&
          //gcpb[iGCPB].dca   > -0.002  && gcpb[iGCPB].dca   < 0.002  &&
          //fabs(gcpb[iGCPB].vert.z) < 100. &&
          gcpb[iGCPB].r_0       >    0     &&
          gcpb[iGCPB].npd_track >    3     &&
          gcpb[iGCPB].p_tot     >    1.e-6 &&
          gcpb[iGCPB].vtl       >    1.e-6 &&
          gcpb[iGCPB].vert.z    > -110.0   && gcpb[iGCPB].vert.z < 110.0 &&
          gcpb[iGCPB].edist     >   -5.0   && gcpb[iGCPB].edist  <  10.0 &&
          gcpb[iGCPB].sdist     >   -5.0   && gcpb[iGCPB].sdist  <   5.0 &&
          gcpb[iGCPB].x2        >    0.3   && gcpb[iGCPB].x2     <   3.0
         )
      {
        rtpc[nRTPC].gcpb   = iGCPB; // pointer to GCPB bank
        rtpc[nRTPC].poverq = gcpb[iGCPB].p_tot;
        rtpc[nRTPC].theta  = gcpb[iGCPB].theta;
        rtpc[nRTPC].phi    = gcpb[iGCPB].phi;
        rtpc[nRTPC].vz     = gcpb[iGCPB].vert.z;
        rtpc[nRTPC].dedx   = gcpb[iGCPB].dedx;
        rtpc[nRTPC].bad    = trackv[iGCPB].bad;

        rtpc[nRTPC].dedx2 = 0;
        rtpc[nRTPC].dedxa = 0;
        rtpc[nRTPC].dedxl = 0.;
        rtpc[nRTPC].dedxal = 0.;

        if (trackv[iGCPB].live_length > 1e-3)
        {
            dx = trackv[iGCPB].live_length;
            rtpc[nRTPC].dedx2  = gcpb[iGCPB].q / dx;
//            rtpc[nRTPC].dedxa  = trackv[iGCPB].qtruncA  / dx;
        }
        if (trackv[iGCPB].live_lengthL > 1e-3)
        {
            dx = trackv[iGCPB].live_lengthL;
            rtpc[nRTPC].dedxl  = trackv[iGCPB].qtruncL  / dx;
        }
        if (trackv[iGCPB].live_lengthA > 1e-3)
        {
            dx = trackv[iGCPB].live_lengthA;
            rtpc[nRTPC].dedxa = trackv[iGCPB].qtruncA / dx;
        }
        if (trackv[iGCPB].live_lengthAL > 1e-3)
        {
            dx = trackv[iGCPB].live_lengthAL;
            rtpc[nRTPC].dedxal = trackv[iGCPB].qtruncAL / dx;
        }

        rtpc[nRTPC].rxy   = trackv[iGCPB].Rxy;
        rtpc[nRTPC].slope = trackv[iGCPB].beta;
        rtpc[nRTPC].chisq = trackv[iGCPB].chisq;
        rtpc[nRTPC].dedxs = trackv[iGCPB].qint / vtl;
        
        // get best hypotheses and energy-loss-corrected momenta:
        eg6rtpc_pids( rtpc[nRTPC].poverq,
                      rtpc[nRTPC].dedx2,
                      rtpc[nRTPC].theta,
                      RTPC_pids,
                      RTPC_moms);
        rtpc[nRTPC].id1 = RTPC_pids[0];
        rtpc[nRTPC].id2 = RTPC_pids[1];
        rtpc[nRTPC].id3 = RTPC_pids[2];
        rtpc[nRTPC].id4 = RTPC_pids[3];
        rtpc[nRTPC].id5 = RTPC_pids[4];
        rtpc[nRTPC].p1  = RTPC_moms[0];
        rtpc[nRTPC].p2  = RTPC_moms[1];
        rtpc[nRTPC].p3  = RTPC_moms[2];
        rtpc[nRTPC].p4  = RTPC_moms[3];
        rtpc[nRTPC].p5  = RTPC_moms[4];

        nRTPC++;

      }
    }

    // fill the RTPC bank:
    if (nRTPC>0)
    {
      clasRTPC_t* RTPC = (clasRTPC_t *) 
        makeBank(&bcs_, "RTPC", 0, sizeof (rtpc_t) / 4, nRTPC);

      for (ii=0; ii < nRTPC; ii++)  RTPC->rtpc[ii]=rtpc[ii];
    }

  } // end fill BOS

  return 1;
}
