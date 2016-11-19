#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <math.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <omp.h> 
using namespace std; 

#define KAPPA 1.0 // diffusivity 
typedef double Real; 

int main(int argc, char *argv[]) {
    // command line argument checking 
    // exactly two arguments required 
    if (argc != 3) {
      printf("USAGE: %s <nx> <nthreads>\n", argv[0]);
      return 1;
    }
    // nx >= 3 required for central difference method
    const int nx = atoi(argv[1]);
    if (nx < 3) { 
        printf("ERROR: nx < 3\n"); 
        return 2; 
    }
    // 0 < nthreads <= nx required 
    const int nthreads = atoi(argv[2]); 
    if (nthreads < 1) { 
        printf("ERROR: nthreads < 1\n"); 
        return 3; 
    }  
    else if (nthreads > nx) { 
        printf("ERROR: nthreads > nx\n"); 
        return 4; 
    } 

    // record program start time (for runtime comparison)
    Real tstart = omp_get_wtime();  

    // set some constants 
    const Real dx = M_PI/nx; // grid spacing
    const Real dt = pow(dx,2)/(8*KAPPA); // time step 
    const Real tf = pow(M_PI,2)/(2*KAPPA); // final time 

    // create objects for parallel file I/O later 
    ostringstream ssnx,ssnth,sstid; 
    ssnx << nx ;
    ssnth << nthreads; 
    string fname; 
    FILE *dataFile; 
    
    // iterator variables 
    int i,j; 

    // initialize data storage matrices  
    Real *Tprev = new Real[nx*nx]; 
    Real *T = new Real[nx*nx];  
    Real Tsum = 0; // variable for volume-integrated average
    Real t = 0; // time variable 

    // store the boundary conditions in arrays 
    // T(x,0) = cos^2(x), T(x,pi) = sin^2(x) 
    Real *upperBoundary = new Real[nx]; 
    Real *lowerBoundary = new Real[nx];  
    for (i=0; i<nx ; i++) { 
        upperBoundary[i] = pow(sin(i*M_PI/nx),2); 
        lowerBoundary[i] = pow(cos(i*M_PI/nx),2); 
    } 

    // temporary storage for points to update with  
    Real Tleft,Tright,Tup,Tdown,Tmid;  
    int tid,jleft,jright;

#pragma omp parallel private(i,j,tid,Tleft,Tright,Tup,Tdown,Tmid,jleft,jright,fname,dataFile,sstid) shared(T,Tprev,upperBoundary,lowerBoundary) reduction(+:Tsum) reduction(+:t) num_threads(nthreads)
    {
        tid=omp_get_thread_num(); 

#pragma omp for schedule(static) 
        // set the initial condition of T=0 everywhere
        for (j=0; j<nx ; j++) { 
            for (i=0; i<nx; i++) { 
                Tprev[j*nx+i]=0;  
            } 
        } 

    // step diffusion equation forward until the final time is reached 
        while (t <= tf) { 
            
#pragma omp for schedule(static)
            for (j=0; j<nx ; j++) { 
                
                // determine the index of the left neighbor 
                if (j==0) { 
                    jleft=nx-1; 
                } 
                else { 
                    jleft=j-1; 
                } 

                // determine the index of the right neighbor 
                if (j==(nx-1)) { 
                    jright=0; 
                } 
                else { 
                    jright=j+1; 
                } 

                for (i=0; i<nx ; i++) { 

                    Tmid=Tprev[j*nx + i]; 
                    Tleft=Tprev[jleft*nx+i]; 
                    Tright=Tprev[jright*nx+i];  
                    
                    // determine value of neighbor below
                    if (i==0) { 
                        Tdown=lowerBoundary[j]; 
                    } 
                    else { 
                        Tdown=Tprev[j*nx + i-1];  
                    } 

                    // determine value of neighbor above 
                    if (i==(nx-1)) { 
                        Tup=upperBoundary[j]; 
                    } 
                    else { 
                        Tup=Tprev[j*nx + i+1];
                    } 
                    
                    // update the i,jth point with centered finited difference in space and forward Euler in time    
                    T[j*nx+i] = Tmid + KAPPA*dt*(Tleft+Tright+Tdown+Tup-4*Tmid)/pow(dx,2); 
                } 
            } 

#pragma omp for schedule(static) 
            // update Tprev and time 
            for (j=0; j<nx ; j++) {
               for (i=0; i<nx; i++) { 
                  Tprev[j*nx+i]=T[j*nx+i]; 
               } 
            }  
            t+=dt; 
        }
    
        // construct file name to write data to
        // each processor writes to its own file 
        sstid << tid; 
        fname = "out_omp." + ssnx.str() + "." + ssnth.str() + "." + sstid.str();
        dataFile = fopen(fname.c_str(),"w"); 

        // compute the final average temperature and write out to console
        // write the final T(x,y) to a file
#pragma omp for schedule(static)
        for (j=0; j<nx; j++) { 
            for (i=0; i<nx; i++) { 
                Tmid=T[j*nx+i];  
                fprintf(dataFile,"%15.8f ",Tmid); 
                Tsum+=Tmid; 
            } 
            fprintf(dataFile,"\n"); 
        } 
        fclose(dataFile);
    } // end parallel section 
    t=t/nthreads; 

    // concatenate all data into single file
    fname = "out_omp." + ssnx.str() + "." + ssnth.str(); 
    char cmd[256]; 
    sprintf(cmd,"a=$(ls %s.*[0-9] | sort -n -t . -k 4); cat $a > %s.all",fname.c_str(),fname.c_str()); 
    system(cmd); 
    // clean up temporary files
    sprintf(cmd,"rm %s.*[0-9]",fname.c_str()); 
    system(cmd); 

    // divide the sum by volume 
    Tsum=Tsum/pow(nx,2); 
   
    // deallocate memory 
    delete [] Tprev; 
    delete [] T; 
    delete [] upperBoundary; 
    delete [] lowerBoundary; 
   
    // determine the program runtime in seconds 
    Real trun=omp_get_wtime()-tstart; 

    // write out summary data to the final data file
    fname=fname + ".all"; 
    dataFile=fopen(fname.c_str(),"a"); 
    fprintf(dataFile,"nx = %d\n",nx); 
    fprintf(dataFile,"trun = %15.8f\n",trun); 
    fprintf(dataFile,"Tavg = %15.8f\n",Tsum); 
    fclose(dataFile); 

    // report successful run execution to console and return
    printf("Average temperature at t=%1.8f: %1.8f\n",t,Tsum); 
    printf("Runtime for nx=%d, nprocs=%d: %1.8f seconds\n",nx,nthreads,trun); 

    return 0; 
} 

