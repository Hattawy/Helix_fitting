#include "rtpc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Number of 8x2 readout boards:
#define __NBOARDS__ 200
// Number of pads per readout board (8x2):
#define __NPADSPERBRD__ 16
// number of noise hits required for rejection:
#define __NNOISEHITS__ 6
// maximum TDC possibly considered noise:
#define __NOISETDCHI__ 32
// number of samples on the noise curve:
#define __NNOISETICS__ (32-12+1)
// number of fired pads on one board above which we apply noise correction:
#define __NPADHOTBOARD__ 10


// count the number of pads that triggered in each board:
void CountBoardPads(int nbrdpads[__NBOARDS__])
{
    int brdpads[__NBOARDS__][__NPADSPERBRD__]; // list of pads in each board
    int ihit,brd,newpad,ii;

    // zero the counters:
    for (ihit=0;ihit<__NBOARDS__;ihit++)
        nbrdpads[ihit]=0;

    for (ihit=0;ihit<hh_num_hits;ihit++)
    {
        if (hh_hitlist[ihit].qraw<=0) continue;
        if (hh_hitlist[ihit].pad<0 || hh_hitlist[ihit].pad>=NUM_PADS) continue;

        if (!IsPadInNoisyBoard(hh_hitlist[ihit].pad)) continue;

        // get board #
        brd=pad2board(hh_hitlist[ihit].pad);

        // check if this pad is already counted:
        newpad=1;
        for (ii=0; ii<nbrdpads[brd]; ii++)
        {
            if (hh_hitlist[ihit].pad == brdpads[brd][ii])
            {
                newpad=0;
                break;
            }
        }

        // add pad to list & increment number of pads in board:
        if (newpad)
            brdpads[brd][nbrdpads[brd]++] = hh_hitlist[ihit].pad;
    }
}
int GetBoardNeighbors(const int board,int ihits[HH_MAX_NUM_HITS])
{
    int nneigh=0;
    if (board<0 || board>=__NBOARDS__) return nneigh;
    
    int row,col,ihit;
    const int row0 = 2*(board/5);
    const int col0 = 8*(board%5);
    const int neighcol[2] = {col0-1,col0+8};
    const int neighrow[2] = {row0-1,row0+2};

    for (ihit=0; ihit<hh_num_hits; ihit++)
    {
        if (hh_hitlist[ihit].qraw<2) continue;

        pad2rowcol(hh_hitlist[ihit].pad,&row,&col);

        if (row==neighrow[0] &&
                col>=neighcol[0] && col<=neighcol[1])
        {
            ihits[nneigh++] = ihit;
        }
        else if (row==neighrow[1] &&
                col>=neighcol[0] && col<=neighcol[1])
        {
            ihits[nneigh++] = ihit;
        }
        else if (col==neighcol[0] &&
                row>=neighrow[0] && row<=neighrow[1])
        {
            ihits[nneigh++] = ihit;
        }
        else if (col==neighcol[1] &&
                row>=neighrow[0] && row<=neighrow[1])
        {
            ihits[nneigh++] = ihit;
        }
    }
    return nneigh;
}
int ApplyHotBoardPedestal(const int board)
{
    // 

    int ii,jj;

    // get indices of hits neighboring the board:
    int neighhits[HH_MAX_NUM_HITS];
    const int nneigh=GetBoardNeighbors(board,neighhits);

    // get indices of hits that are in the board:
    int nboardhits=0;
    int boardhits[HH_MAX_NUM_HITS];
    for (ii=0; ii<hh_num_hits; ii++)
    {
        if (pad2board(hh_hitlist[ii].pad) == board)
        {
            boardhits[nboardhits++] = ii;
        }
    }

    // calculate pedestal from hits on board without off-board neighbors: 
    int nped=0;
    float qq,pedestal=0;
    int row,col,hasneighbor;
    int neighrow,neighcol;
    for (ii=0; ii<nboardhits; ii++) // hits on hot board
    {
        pad2rowcol(hh_hitlist[boardhits[ii]].pad,&row,&col);

        hasneighbor=0;
        for (jj=0; jj<nneigh; jj++) // hits on board neighbors
        {
            pad2rowcol(hh_hitlist[neighhits[jj]].pad,&neighrow,&neighcol);

            if (abs(neighcol-col) <=1)
            {
                hasneighbor=1;
                break;
            }
        }
        // no off-board neighbors, include in pedestal calculation:
        if (!hasneighbor)
        {
            nped++;
            pedestal += hh_hitlist[boardhits[ii]].qraw;
        }
    }

    // insufficient hit without neighbors, no pedestal correction:
    if (nped>2) pedestal /= nped;
    else return 0;

    // apply pedestal to all hits on the board (RESET hh_hitlist.q !!!):
    for (ii=0; ii<nboardhits; ii++)
    {
        qq = hh_hitlist[boardhits[ii]].qraw;
        qq += altro_offset - pedestal;
        qq /= corr_fac[hh_hitlist[boardhits[ii]].pad];
        if (qq < altro_offset+5) 
        {
            qq=0;
            hh_hitlist[boardhits[ii]].status |= HBADPAD;
        }
        hh_hitlist[boardhits[ii]].q = qq;
    }
    return nboardhits;
}
int RejectBoard(const int board)
{
    int ihit,nrej=0;
    for (ihit=0; ihit<hh_num_hits; ihit++)
    {
        if (hh_hitlist[ihit].pad<0 || hh_hitlist[ihit].pad>=NUM_PADS) 
            continue;
        if (pad2board(hh_hitlist[ihit].pad) == board)
        {
            hh_hitlist[ihit].status |= HBADPAD;
            nrej++;
        }
    }
    return nrej;
}
int RejectBusyBoards()
{
     int nrej=0,ihit,nbrdpads[__NBOARDS__];
     CountBoardPads(nbrdpads);
     for (ihit=0; ihit<hh_num_hits; ihit++)
     {
        if (hh_hitlist[ihit].pad<0 || hh_hitlist[ihit].pad>=NUM_PADS) continue;
        if (nbrdpads[pad2board(hh_hitlist[ihit].pad)] > __NPADHOTBOARD__)
        {
            hh_hitlist[ihit].status |= HBADPAD;
            nrej++;
        }
     }
     return nrej;
}
int BusyBoards()
{
     int brd,nrej=0;
     int nbrdpads[__NBOARDS__];

     // count number of pads triggered on each board:
     CountBoardPads(nbrdpads);

     for (brd=0; brd<__NBOARDS__; brd++)
     {
         if (nbrdpads[brd] <= __NPADHOTBOARD__) continue;

         // if pedestal correction is not possible, reject all board's hits:
         if (!ApplyHotBoardPedestal(brd))
             nrej += RejectBoard(brd); 
     }
     return nrej;
}

int RejectNoise()
{
    int nnoise=0;
    int ipad,ihit,nhit;
    int tdcs[NAL_SAMP],adcs[NAL_SAMP],hits[NAL_SAMP];

    for (ipad=0; ipad<NUM_PADS; ipad++)
    {
        // don't bother if pad is marked bad or not noisy:
        if (dead_pad(ipad)) continue;
        if (!IsPadNoisy(ipad)) continue;
//        if (!IsPadInNoisyBoard(ipad)) continue;

        // load arrays of hits for one pad:
        nhit=0;
        for (ihit=0; ihit<hh_num_hits; ihit++)
        {
            if (hh_hitlist[ihit].pad != ipad) continue;

            // don't bother if TDC is above noise region:
            if (hh_hitlist[ihit].t > __NOISETDCHI__) continue;

            hits[nhit] = ihit;
            tdcs[nhit] = hh_hitlist[ihit].t;
            adcs[nhit] = hh_hitlist[ihit].qraw;
            nhit++;
           
            // this should never happen, but:
            if (nhit >= NAL_SAMP) break;
        }

        // count number of noise hits:
        if (IsNoise(nhit,tdcs,adcs)<__NNOISEHITS__) continue;
       
        // mark hits as noise, bad:
        for (ihit=0; ihit<nhit; ihit++)
        {
            if (tdcs[ihit] > __NOISETDCHI__) continue;
            hh_hitlist[hits[ihit]].status |= HBADPAD;
            nnoise++;
        }
    }
    return nnoise;
}
int IsOnNoiseCurve2(const int tdc,const int adc)
{
    // check whether hit lies on noise curve.
 
//    // ORIGINAL:
////  static int tdcs[__NNOISETICS__]={12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
//    static int adcs[__NNOISETICS__]={ 0, 0, 0,30,67,65,43,40,25, 0, 0,30,67,65,45,42,28,26, 0, 0, 0};
//    static int eeee[__NNOISETICS__]={15,20,25,15,17,10,18,21,15,20,25,15,17,10,18,20,15,14, 0, 0, 0};

//    // INCREASE WIDTHS on 15,17,23,25:
////  static int tdcs[__NNOISETICS__]={12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
//    static int adcs[__NNOISETICS__]={ 0, 0, 0,30,67,65,43,40,25, 0, 0,30,67,65,45,42,28,26, 0, 0, 0};
//    static int eeee[__NNOISETICS__]={15,20,25,23,23,23,18,21,15,20,25,23,23,23,18,20,15,14, 0, 0, 0};

    // SHIFT 15,23:
//  static int tdcs[__NNOISETICS__]={12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    static int adcs[__NNOISETICS__]={ 0, 0, 0,50,67,65,43,40,25, 0, 0,50,67,65,45,42,28,26, 0, 0, 0};
    static int eeup[__NNOISETICS__]={15,20,25,23,23,23,18,21,15,20,25,23,23,23,18,20,15,14, 0, 0, 0};
    static int eedn[__NNOISETICS__]={-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0};
    
    // adcs is the average ADC value of noise
    // eeup is the ADC spread of noise
    // eedn is the lower limit.  If negative, ignored and using spread instead.

    const int tdcbin=tdc-12;
    if (tdcbin<0 || tdcbin>=__NNOISETICS__) return 0;

    // ADC way too big for pure noise:
    if (adc > adcs[tdcbin] + 2*eeup[tdcbin]) return -1;
    // ADC too small for pure noise:
    if (eedn[tdcbin]>=0 && adc < adcs[tdcbin] - 1.5*eeup[tdcbin]) return -1; 

    // ADC too big for pure noise:
    if (adc > adcs[tdcbin]+eeup[tdcbin]) return 0;

    // ADC too small for pure noise:
    if (eedn[tdcbin] >= 0) 
    {
        if (adc < eedn[tdcbin]) return 0;
    }
    else if (adc < adcs[tdcbin]-eeup[tdcbin]) return 0;

    return 1;

//    return fabs(adc-adcs[tdcbin])<eeee[tdcbin];
}
int IsNoise(int nhit,const int tdc[NAL_SAMP],const int adc[NAL_SAMP])
{
    // count number of hits that fall on the noise curve

    if (nhit>NAL_SAMP) nhit=NAL_SAMP;
    int ii,nnoise=0,nnoise1=0,nnoise2=0,nbig=0;
    for (ii=0; ii<nhit; ii++)
    {
        // count big hits:
        if (IsOnNoiseCurve2(tdc[ii]  ,adc[ii])<0) nbig++;
        
        // count noise hits:
        // allow +/- 1 shift in TDC
        if (IsOnNoiseCurve2(tdc[ii]  ,adc[ii])>0) nnoise++;
        if (IsOnNoiseCurve2(tdc[ii]-1,adc[ii])>0) nnoise1++;
        if (IsOnNoiseCurve2(tdc[ii]+1,adc[ii])>0) nnoise2++;
    }
    if (nnoise1 > nnoise) nnoise = nnoise1;
    if (nnoise2 > nnoise) nnoise = nnoise2;
    
    return nbig>1 ? -nnoise : nnoise;
}

#define __NNOISYPADS__ 581
int IsPadNoisy(const int pad)
{
  if (pad<0 || pad>=3200) return 0;

  static const int noisypads[__NNOISYPADS__]=
  { 150,151,152,153,154,155,156,157,158,159,160,161,
    162,163,164,165,166,167,168,169,170,171,172,173,
    174,248,250,251,252,253,274,275,276,277,278,279,
    281,282,283,284,285,286,287,400,401,402,403,404,
    405,406,407,408,409,410,411,412,413,414,415,416,
    417,418,419,420,421,422,423,424,425,426,427,428,
    429,430,431,433,434,435,436,437,438,439,445,446,
    447,487,532,533,534,535,536,537,538,539,540,541,
    542,543,656,657,658,659,660,661,662,663,664,665,
    666,667,668,669,670,671,
    1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,
    1050,1051,1052,1053,1054,1055,1168,1169,1170,1171,
    1172,1173,1174,1175,1176,1177,1178,1179,1180,1181,
    1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,
    1192,1193,1194,1195,1196,1197,1198,1199,1297,1297,
    1299,1300,1301,1302,1303,1304,1305,1306,1307,1308,
    1309,1310,1311,1314,1315,1316,1317,1321,1323,1324,
    1325,1365,1366,1367,1374,1375,1428,1429,1430,1431,
    1432,1433,1434,1435,1436,1437,1438,1439,1744,1745,
    1746,1747,1749,1751,1752,1753,1755,1758,1760,1761,
    1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,
    1772,1773,1774,1775,1878,1879,1880,1881,1882,1883,
    1884,1885,1886,1887,1889,1890,1891,1892,1893,1894,
    1895,1896,1897,1898,1899,1900,1901,1902,1903,1907,
    1908,1909,1910,1911,1912,1917,1918,1919,1920,1921,
    1922,1923,1924,1925,1926,1927,1928,1970,1971,1972,
    1973,1974,1975,1976,1977,1981,1982,1983,1984,1985,
    1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,
    1996,1997,1998,1999,2032,2033,2034,2035,2036,2037,
    2038,2039,2045,2046,2047,2048,2049,2050,2051,2052,
    2053,2054,2055,2061,2062,2063,2128,2129,2130,2131,
    2131,2132,2134,2135,2136,2137,2138,2139,2140,2141,
    2142,2143,2144,2145,2146,2147,2148,2149,2150,2151,
    2152,2153,2154,2155,2156,2157,2158,2159,2232,2245,
    2246,2247,2248,2249,2250,2253,2254,2255,2256,2257,
    2258,2259,2260,2261,2262,2263,2264,2265,2266,2267,
    2268,2269,2270,2271,2272,2273,2274,2275,2276,2277,
    2278,2279,2280,2281,2282,2283,2284,2285,2286,2287,
    2368,2369,2370,2371,2372,2373,2374,2375,2376,2377,
    2378,2379,2380,2381,2382,2383,2391,2396,2397,2398,
    2399,2401,2402,2403,2404,2405,2406,2407,2408,2409,
    2410,2411,2412,2413,2414,2498,2499,2500,2501,2502,
    2503,2504,2505,2506,2509,2510,2511,2512,2513,2514,
    2515,2516,2517,2518,2519,2520,2521,2522,2523,2524,
    2525,2526,2527,2528,2529,2530,2531,2532,2533,2534,
    2535,2536,2537,2538,2539,2540,2541,2542,2543,2552,
    2640,2641,2642,2643,2644,2645,2646,2647,2648,2649,
    2650,2651,2652,2653,2654,2655,2768,2769,2770,2771,
    2772,2773,2774,2775,2776,2777,2778,2779,2780,2781,
    2782,2783,2784,2785,2786,2787,2788,2789,2790,2791,
    2792,2793,2794,2795,2796,2797,2798,2799,2896,2896,
    2897,2898,2899,2900,2901,2902,2903,2904,2905,2906,
    2907,2908,2909,2910,2911,2912,2913,2914,2915,2916,
    2917,2918,2919,2920,2921,2922,2923,2924,2925,2926,
    2927,3015,3016,3017,3018,3021,3024,3025,3026,3027,
    3028,3029,3030,3031,3032,3033,3034,3035,3036,3037,
    3038,3039,3195,3196,3197,3198,3199};


  int ii,imin,imax;
  static int noisypads2[NUM_PADS]={}; // all zeroes by default
  static int loaded=0;

  // load array of 3200 noisy flags based upon pad list:
  if (!loaded)
  {
      for (ii=0; ii<NUM_PADS; ii++)
      {

          // binary search:
          imin=0;
          imax=__NNOISYPADS__-1;
          while (imax >= imin)
          {
              const int imid=(imin+imax)/2;
              if (noisypads[imid] == ii) 
              {
                  noisypads2[ii]=1;
                  break;
              }
              else if (noisypads[imid]  < pad) imin = imid+1;
              else                             imax = imid-1;
          }
      }
      loaded=1;
  }

  return noisypads2[pad];

/*
  // binary search:
  int imin=0;
  int imax=__NNOISYPADS__-1;
  while (imax >= imin)
  {
      const int imid=(imin+imax)/2;
      if      (noisypads[imid] == pad) return 1;
      else if (noisypads[imid]  < pad) imin = imid+1;
      else                             imax = imid-1;
  }
  return 0;
*/
}
#define __NNOISYBOARDS__ 45
int IsPadInNoisyBoard(const int pad)
{
  if (pad<0 || pad>=NUM_PADS) return 0;

  static const int noisybrds[__NNOISYBOARDS__]=
  {  9, 10, 15, 17, 25, 26, 27, 33, 41, 65, 73,
    74, 81, 82, 85, 89,109,110,117,118,119,120,
   123,124,127,128,133,134,140,141,142,148,149,
   150,156,157,158,165,173,174,181,182,188,189,199 };

  int ii,imin,imax,brd;
  static int noisybrdspads[NUM_PADS]={};
  static int loaded=0;

  // load array of 3200 noisy flags based upon board list:
  if (!loaded)
  {
      for (ii=0; ii<NUM_PADS; ii++)
      {
          brd=pad2board(ii);

          // binary search:
          imin=0;
          imax=__NNOISYBOARDS__-1;
          while (imax >= imin)
          {
              const int imid=(imin+imax)/2;
              if (noisybrds[imid] == brd) 
              {
                  noisybrdspads[ii]=1;
                  break;
              }
              else if (noisybrds[imid]  < brd) imin = imid+1;
              else if (noisybrds[imid]  > brd) imax = imid-1;
          }
      }
      loaded=1;
  }

  return noisybrdspads[pad];

/*
  if (pad<0 || pad>=NUM_PADS) return 0;

  const int brd=pad2board(pad);
 
  // binary search:
  int imin=0;
  int imax=__NNOISYBOARDS__-1;
  while (imax >= imin)
  {
      const int imid=(imin+imax)/2;
      if      (noisybrds[imid] == brd) return 1;
      else if (noisybrds[imid]  < brd) imin = imid+1;
      else if (noisybrds[imid]  > brd) imax = imid-1;
  }
  return 0;
*/
}
