#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <math.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
using namespace std; 

#define KAPPA 1.0 // diffusivity 
typedef double Real; 

int main(int argc, char *argv[]) {
    // command line argument checking 
    // exactly one argument required 
    if (argc != 2) {
        printf("USAGE: %s <nx>\n", argv[0]);
        return 1;
    }
    // nx >= 3 required for central difference method
    const int nx = atoi(argv[1]);
    if (nx < 3) { 
        printf("ERROR: nx < 3\n"); 
        return 2; 
    }

    // set some constants 
    const Real dx = M_PI/nx; // grid spacing
    const Real dt = pow(dx,2)/(8*KAPPA); // time step 
    const Real tf = pow(M_PI,2)/(2*KAPPA); // final time 

    // initialize MPI 
    int rc,ntasks,rank; 
    rc = MPI_Init(&argc,&argv); 
    if (rc != MPI_SUCCESS) { 
        printf("ERROR: MPI_Init failed"); 
        MPI_Abort(MPI_COMM_WORLD,rc); 
        return 3; 
    } 
    MPI_Comm_size(MPI_COMM_WORLD,&ntasks); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // determine size for each thread and offset (for boundary condition) 
    int njper=nx/ntasks;   
    int mod = fmod(njper,1); 
    if (rank < mod) { 
        njper++;
    } 
    int jbeg = rank*njper; 
    if (rank >= mod) { 
        jbeg+= mod; 
    } 
    
    // record start time for timing 
    double tstart = MPI_Wtime(); 

    // nxfull is the full size of each array, accounting for 2 ghost points
    const int njfull = njper + 2; 
    const int ni = nx; 

    // set ranks of nearest neighbors 
    int rankLeft, rankRight; 
    if (rank == 0) { 
        rankLeft=ntasks-1; 
    } 
    else { 
        rankLeft=rank-1; 
    } 
    if (rank == (ntasks-1)) { 
        rankRight=0; 
    } 
    else { 
        rankRight=rank+1;
    } 

    // initialize data storage matrices  
    Real *Tprev = new Real[ni*njfull]; 
    Real *T = new Real[ni*njfull]; 
    Real Tsum = 0; // variable for volume-integrated average
    Real t = 0; // time variable 

    // store the boundary conditions and set the initial conditions  
    Real *upperBoundary = new Real[njper]; 
    Real *lowerBoundary = new Real[njper];  
    int i,j; 
    for (j=0; j < njfull ; j++) { 
        if (j < njper) { 
            // boundary conditions: T(x,0) = cos^2(x), T(x,pi) = sin^2(x)
            upperBoundary[j] = pow(sin((jbeg+j)*M_PI/nx),2); 
            lowerBoundary[j] = pow(cos((jbeg+j)*M_PI/nx),2); 
        } 
        // initial conditions: T(x,y) = 0
        for (i=0; i<ni ; i++) { 
            Tprev[j*ni + i] = 0; 
        } 
    } 
    
    // temporary storage for points to update with  
    Real Tleft,Tright,Tup,Tdown,Tmid;  
    MPI_Status Stat; 
    MPI_Request Req1,Req2; 
// step diffusion equation forward until the final time is reached 
    while (t <= tf) { 

        // update physical grid points
        // (j=0 and j=njfull-1 are ghost points)
        for (j=1; j<njfull-1 ; j++) {

            for (i=0; i<ni; i++) { 
                Tmid=Tprev[j*ni + i]; 
                Tleft=Tprev[(j-1)*ni + i]; 
                Tright=Tprev[(j+1)*ni + i]; 
                
                // determine value of neighbor below
                // edge cases have j-1 instead of j because loop starts at j=1
                if (i==0) { 
                    Tdown=lowerBoundary[j-1]; 
                } 
                else {
                    Tdown=Tprev[j*ni + i-1];  
                } 

                // determine value of neighbor above 
                if (i==(nx-1)) { 
                    Tup=upperBoundary[j-1]; 
                } 
                else { 
                    Tup=Tprev[j*ni + i+1]; 
                } 
                
                // update the i,jth point with centered finited difference in space and forward Euler in time    
                T[j*ni+i]=Tmid + KAPPA*dt*(Tleft+Tright+Tdown+Tup-4*Tmid)/pow(dx,2); 
            } 
        } 

        // communicate ghost points 
        // send ghost points to the left and right neighbors
        MPI_Isend(&(T[ni]),ni,MPI_DOUBLE,rankLeft,0,MPI_COMM_WORLD,&Req1); 
        MPI_Isend(&(T[(njfull-2)*ni]),ni,MPI_DOUBLE,rankRight,2,MPI_COMM_WORLD,&Req2); 
       
        // wait for ghost points from the left and right neighbors  
        MPI_Recv(&(T[0]),ni,MPI_DOUBLE,rankLeft,2,MPI_COMM_WORLD,&Stat);  
        MPI_Recv(&(T[(njfull-1)*ni]),ni,MPI_DOUBLE,rankRight,0,MPI_COMM_WORLD,&Stat);  

        // update Tprev and time 
        for (j=0; j<njfull ; j++) { 
            for (i=0; i<ni ; i++) { 
                Tprev[j*ni + i]=T[j*ni+i]; 
            } 
        }
        t+=dt; 
    } // end time evolution 

    // prepare file for writing 
    // each processor writes to its own file 
    ostringstream ssnx,ssnth,sstid; 
    ssnx << nx ;
    ssnth << ntasks; 
    sstid << rank; 
    string fname = "out_mpi." + ssnx.str() + "." + ssnth.str() + "." + sstid.str();
    FILE *dataFile = fopen(fname.c_str(),"w"); 

    // compute the final average temperature and write out to console
    // write the final T(x,y) to a file
    for (j=1; j<njfull-1; j++) { 
        for (i=0; i<ni; i++) { 
            Tmid = T[j*ni + i]; 
            fprintf(dataFile,"%15.8f ",Tmid); 
            Tsum+=Tmid; 
        } 
        fprintf(dataFile,"\n"); 
    } 
    fclose(dataFile);
    
    // add together sums from all processors 
    Real Tsumall=0;
    MPI_Allreduce(&Tsum,&Tsumall,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
    
    // deallocate memory 
    delete [] Tprev; 
    delete [] T; 
    delete [] upperBoundary; 
    delete [] lowerBoundary; 
   
    // master process concatenates all data files together 
    // and writes data to console 
    if (rank == 0) { 
        Tsumall=Tsumall/pow(nx,2); 

        // concatenate all data into single file
        fname = "out_mpi." + ssnx.str() + "." + ssnth.str(); 
        char cmd[256]; 
        sprintf(cmd,"a=$(ls %s.*[0-9] | sort -n -t . -k 4); cat $a > %s.all",fname.c_str(),fname.c_str()); 
        system(cmd);
        // clean up temporary files
        sprintf(cmd,"rm %s.*[0-9]",fname.c_str()); 
        system(cmd); 
        // runtime elapsed 
        double trun = MPI_Wtime()-tstart; 

        // write out summary data to the final data file
        fname=fname + ".all"; 
        dataFile=fopen(fname.c_str(),"a"); 
        fprintf(dataFile,"nx = %d\n",nx); 
        fprintf(dataFile,"trun = %15.8f\n",trun); 
        fprintf(dataFile,"Tavg = %15.8f\n",Tsumall); 
        fclose(dataFile); 

        // report successful run execution to console and return
        printf("Average temperature at t=%1.8f: %1.8f\n",t,Tsumall); 
        printf("Runtime for nx=%d, nprocs=%d: %1.8f seconds\n",nx,ntasks,trun); 
    } 
   
    // Finalize MPI
    MPI_Finalize(); 

    return 0; 
} 

