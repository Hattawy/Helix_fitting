// EG6 RTPC Energy Loss Correction
// by Michael Oliver
//
// Correction for energy lost between vertex and drift region.
// Data tables are calculated from simulation.
// Correction function does 2-D interpolation in p/theta.
//
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define BUFSIZE 256
#define ERR_MISSING_FILE 1
#define ERR_BAD_ASCII_FILE 2
#define ERR_OUTPUT 0

// planar interpolation algorithm:
#define PLANAR(x1,y1,z1,x2,y2,z2,x3,y3,z3,x,y) (DET(x1,y1,z1,x2,y2,z2,x3,y3,z3) - x * DET(1,y1,z1,1,y2,z2,1,y3,z3) - y * DET(x1,1,z1,x2,1,z2,x3,1,z3) ) / DET(x1,y1,1,x2,y2,1,x3,y3,1)
#define DET(a1,b1,c1,a2,b2,c2,a3,b3,c3) (a1 * (b2*c3 - c2*b3) - a2 * (b1*c3 - c1*b3) + a3 * (b1*c2 - c1*b2))

// number of recoil types (p,d,T,3He,4He):
#define NPARTS 5

// BEWARE: these are hardcoded to match the input data table. ///////////
// step size:
#define P_STEP 5
#define T_STEP 1
// array lengths:
#define P_LENGTH 300
#define T_LENGTH 181 
// momentum cutoff, no energy loss above this (p,d,T,3He,4He):
#define MAX_P_ARRAY {805,945,1080,1255,1410}
/////////////////////////////////////////////////////////////////////////


/**
 * loads the data tables -- only called one time
 * first try $CLAS_PACK/eg6eloss_data/, then ./eg6eloss_data/ 
 */
int eg6rtpc_elossload (float data_array[NPARTS][T_LENGTH][P_LENGTH]) 
{
    const char *stubs[NPARTS] = {"proton.txt","deuteron.txt","triton.txt","3He.txt","4He.txt"};
    const char *prefix = "eg6eloss_data/eg6eloss";
    FILE *file; 
    char buf[BUFSIZE],cc,filename[1000];
    int num_lines,ii,pid_n,ix,iy;
    double p_i,t_i,p_f;

    for (pid_n = 0; pid_n < NPARTS; pid_n++)
    {
        // try current directory first:
        sprintf(filename,"./%s_%s",prefix,stubs[pid_n]);
        if (NULL == (file = fopen (filename, "r")))
        {
            // then try CLAS_PARMS: 
            sprintf(filename,"%s/%s_%s",getenv("CLAS_PARMS"),prefix,stubs[pid_n]);
            if (NULL == (file = fopen (filename, "r")))
            {
                fprintf(stderr,"EG6RTPCELOSS Error: missing file \"%s\"\n",filename);
                return ERR_MISSING_FILE;
            }
        }
        fprintf(stdout,"EG6RTPCELOSS reading data: %s\n",filename);

        // get number of lines in file:
        num_lines = 0;
        while ((cc = getc(file)) != EOF) {
            if (cc == '\n') num_lines++;
        }
        rewind(file);
        
        // read the file:
        ii = 0;
        while ((fgets(buf,BUFSIZE,file)) != NULL)
        {    
            if (ii > num_lines ) 
            {
                fprintf(stderr,"EG6RTPCELOSS Error: bad ascii file \"%s\"\n",filename);
                return ERR_BAD_ASCII_FILE;
            }
           
            p_f = atof( (char *) strtok(buf," \t")    );
            t_i = atof( (char *) strtok(NULL," \t")   );
            p_i = atof( (char *) strtok(NULL," \t\n") );

            if (isnan(p_i) || isnan(p_f) || isnan(t_i)) 
            {
                fclose(file);
                fprintf(stderr,"EG6RTPCELOSS Error: bad ascii file \"%s\"\n",filename);
                return ERR_BAD_ASCII_FILE;
            }
           
            // get array indices:
            iy = (int) (p_f / P_STEP);
            ix = (int) (t_i / T_STEP);
      
            if (iy<0 || ix<0 || iy>=P_LENGTH || ix>=T_LENGTH)
            {
                fprintf(stderr,"EG6RTPCELOSS Error: ascii file data out of array bounds.");
                return ERR_BAD_ASCII_FILE;
            }

            data_array[pid_n][ix][iy] = p_i;

            ii++;
        }
        
        fclose(file);

        if (ii!= num_lines)
        {
            fprintf(stderr,"EG6RTPCELOSS Error: mismatch in ascii file \"%s\"\n",filename);
            return ERR_BAD_ASCII_FILE;
        }

    }
    return 0;
}


/**
 * Given particle type and theta and real momentum in the drift
 * region (not p/q), return real vertex-momentum by interpolating
 * table from simulation.  Load the table if necessary.
 *
 * Input pid is 2212,45,46,49,47 for p,d,T,3He,4He.
 * Input theta and momentum are in radians and MeV.
 */
float eg6rtpc_eloss (int pid, float t_i, float p_f)
{
    static float data_array[NPARTS][T_LENGTH][P_LENGTH]={};
    static char data_loaded = 0;
    static char data_good = 0;
    const float max_momenta[NPARTS] = MAX_P_ARRAY;
    double c[5];        // closest bins' vertex-momenta
    double t_k,t_k_1;   // theta kk, theta kk+1
    double pf_m,pf_m_1; // drift-momentum in bins mm, mm+1, and mm+2
    int pid_n;          // 0,1,2,3,4 = p,d,T,3He,4He
    int mm,kk;          // p and theta bin numbers
    int i_p,i_t;        // closest bin flags
    double result;
    
    // load data:
    if (!data_loaded)
    {
        data_loaded = 1;
        if(!eg6rtpc_elossload(data_array)) data_good = 1;
        else fprintf(stderr, "EG6RTPCELOSS Error:  Bad Data Read"
                             " ->  Disabling Correction !!!!!!!!!!!\n");
    }
    if (!data_good) return ERR_OUTPUT;

    // convert input GEANT3 pid to array index:
    switch (pid) 
    {
        case 2212:     //proton
            pid_n = 0;
            break;
        case 45:       //deuteron
            pid_n = 1;
            break;
        case 46:       //triton
            pid_n = 2;
            break;
        case 49:       //3He
            pid_n = 3; 
            break;
        case 47:       //4He
            pid_n = 4;
            break;
        default:
            fprintf(stderr,"EG6RTPCELOSS Error: unkown pid %d\n",pid);
            return ERR_OUTPUT;
    }

    // data files are in degrees:
    t_i *= 57.29578;

    // check for kinematics outside valid range:
    if (fabs(t_i-90) >= 78 || p_f < 10) return ERR_OUTPUT;

    // no energy loss correction above maximum momentum:
    if (p_f + 10 >= max_momenta[pid_n]) return p_f;
    
    // correction defined only for 90-180, but it's symmetric about 90:
    if (t_i < 90) t_i = 180 - t_i;
    
    // theta and drift-momentum bin numbers:
    kk = (int) (t_i/T_STEP);
    mm = (int) (p_f/P_STEP);

    // theta and drift-momentum for this bin:
    t_k  = T_STEP * kk;
    pf_m = P_STEP * mm;
   
    // vertex-momenta of the four closest points in the table:
    c[0] = data_array[pid_n][kk][mm];
    c[1] = data_array[pid_n][kk][mm+1];
    c[2] = data_array[pid_n][kk+1][mm];
    c[3] = data_array[pid_n][kk+1][mm+1];

    //fprintf(stderr,"t/p=%6.2f/%6.2f   ",t_i,p_f);
    //fprintf(stderr,"k/m=%d/%d <==> t/p=%4.0f/%4.0f   ",kk,mm,t_k,pf_m);
    //fprintf(stderr,"0/1/2/3=%6.2f/%6.2f/%6.2f/%6.2f   ",*c[0],*c[1],*c[2],*c[3]);

    // check for hole in data (fiducial cut will remove this case):
    if (!(c[0]*c[1]*c[2]*c[3])) {
        // moving in positive p until no hole:
        mm=0;
        do {
            mm++;
            if (mm > 100) return ERR_OUTPUT;
            c[0] = data_array[pid_n][kk][mm];
            c[1] = data_array[pid_n][kk][mm+1];
            c[2] = data_array[pid_n][kk+1][mm];
            c[3] = data_array[pid_n][kk+1][mm+1];
            c[4] = data_array[pid_n][kk][mm+2];
        }
        while (!(c[0]*c[1]*c[2]*c[3]*c[4]));
	return c[0];
    }

    // theta & drift-momenta of the closest 4 points:
    t_k    = (double)(T_STEP * kk);
    t_k_1  = (double)(T_STEP * (kk + 1));
    pf_m   = (double)(P_STEP * mm);
    pf_m_1 = (double)(P_STEP * (mm + 1));

    // flags to find the 3 closest points:
    //i_t = 0 => t_i closest to t_k
    //i_t = 1 => t_i closest to t_k_1
    //i_p = 0 => p_f closest to pf_m
    //i_p = 1 => p_f closest to pf_m_1
    i_t = ((t_i - t_k) / T_STEP) < 0.5;
    i_p = ((p_f - pf_m) / P_STEP) < 0.5;

    // do interpolation using the closest 3 points:
    switch (i_p + 2*i_t) 
    {
        case 0:
            result = PLANAR(pf_m_1, t_k_1, c[3],
                            pf_m_1, t_k,   c[1],
                            pf_m,   t_k_1, c[2],
                            (double)p_f, (double)t_i);
            break;
        case 1:
            result = PLANAR(pf_m, t_k, c[0],
                            pf_m_1, t_k_1, c[3],
                            pf_m,   t_k_1, c[2],
                            (double)p_f, (double)t_i);
            break;
        case 2:
            result = PLANAR(pf_m, t_k, c[0],
                            pf_m_1, t_k,   c[1],
                            pf_m_1, t_k_1, c[3],
                            (double)p_f, (double)t_i);
            break;
        default:  //3
            result = PLANAR(pf_m, t_k, c[0],
                            pf_m_1, t_k,   c[1],
                            pf_m,   t_k_1, c[2],
                            (double)p_f, (double)t_i);
    }
    
    return isnan(result) ? ERR_OUTPUT : result;
}



/*
// to compile standalone for testing:
int main (int argc, char **argv)
{ 
    double result;
    double t_i, p_f;
    int pid; //particle id
    char buf[BUFSIZE];
    char *str;
    const int pids[5]={2212,45,46,49,47};

    printf("Input: <particle id> <initial theta> <final momentum>\n");
    while (1)
    {    
        if (NULL==fgets(buf,BUFSIZE,stdin)) break; 
        str = (char *) strtok(buf," ");
        if (str[0] == '\n') 
        {
            // infinite loop with random input:
            srand(9439237);
            while (1)
            {
                pid=pids[rand()%5];
                t_i=180*(float)rand()/RAND_MAX;
                p_f=2000*(float)rand()/RAND_MAX;
                result = eg6rtpc_eloss(pid,t_i,p_f);
                printf("%d %f %f %f\n",pid,t_i,p_f,result);
            }
        }
        if (buf[0] == 'q') exit(0);
        
        pid = atoi(str);
        
        str = (char *) strtok(NULL," ");
        t_i = atof(str);
        str =  (char *) strtok(NULL," ");
        p_f = atof(str);

        result = eg6rtpc_eloss(pid,t_i,p_f);
        printf("%d %f %f %f\n",pid,t_i,p_f,result);
    }
    exit(0);
}
*/

