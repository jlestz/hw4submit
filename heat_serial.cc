#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
#include <math.h> 
#include <iostream> 
#include <fstream> 
#include <sstream> 
#include <time.h> 
using namespace std; 

#define KAPPA 1.0 // diffusivity 
typedef double Real; 

int main(int argc, char *argv[]) {
    // command line argument checking 
    // nx >= 3 required for central difference method
    if (argc != 2) {
      printf("USAGE: %s <nx>\n", argv[0]);
      return 1;
    }
    const int nx = atoi(argv[1]);
    if (nx < 3) { 
        printf("ERROR: nx < 3\n"); 
        return 2; 
    }

    // record program start time (for runtime comparison)
    Real tstart = clock();  

    // set some constants 
    const Real dx = M_PI/nx; // grid spacing
    const Real dt = pow(dx,2)/(8*KAPPA); // time step 
    const Real tf = pow(M_PI,2)/(2*KAPPA); // final time 

    // initialize data storage matrices  
    Real *Tprev = new Real[nx*nx]; 
    Real *T = new Real[nx*nx]; 
    Real *upperBoundary = new Real[nx]; 
    Real *lowerBoundary = new Real[nx];  

    // set the initial condition and boundary conditions 
    int i,j; 
    for (j=0; j<nx ; j++) { 
       for (i=0; i<nx ; i++) {
          // initial condition T(x,y) = 0
          Tprev[j*nx + i] = 0; 
       }
       // boundary conditions T(x,0) = cos^2(x), T(x,pi) = sin^2(x) 
       upperBoundary[j] = pow(sin(j*M_PI/nx),2); 
       lowerBoundary[j] = pow(cos(j*M_PI/nx),2); 
    }

    // temporary storage for points to update with  
    Real Tleft,Tright,Tup,Tdown,Tmid;  
    int jleft,jright; 

    // step diffusion equation forward until the final time is reached 
    Real t = 0; 
    while (t <= tf) { 
        
        for (j=0; j<nx ; j++) { 

            // determine index of the left neighbor 
            if (j==0) { 
                jleft=nx-1; 
            } 
            else { 
                jleft=j-1; 
            } 

            // determine index of the right neighbor 
            if (j==(nx-1)) { 
                jright=0; 
            } 
            else { 
                jright=j+1; 
            } 

            for (i=0; i<nx ; i++) { 

                Tmid=Tprev[j*nx + i]; 
                Tleft=Tprev[jleft*nx + i]; 
                Tright=Tprev[jright*nx + i]; 

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

        // update Tprev and time 
        for (j=0; j<nx; j++) { 
            for (i=0; i<nx; i++) { 
                Tprev[j*nx+i]=T[j*nx+i]; 
            } 
        } 
        t+=dt; 
    } 
    
    // construct file name to write data to
    ostringstream ss; 
    ss << nx ;
    string fname = "out_serial." + ss.str(); 
    FILE *dataFile = fopen(fname.c_str(),"w"); 

    // compute the final average temperature and write out to console
    // write the final T(x,y) to a file
    Real sum = 0; 
    for (j=0; j<nx; j++) { 
        for (i=0; i<nx; i++) { 
            Tmid=T[j*nx + i]; 
            fprintf(dataFile,"%15.8f ",Tmid); 
            sum+=Tmid; 
        } 
        fprintf(dataFile,"\n"); 
    } 
    sum=sum/pow(nx,2); 
   
    // deallocate memory 
    delete [] Tprev; 
    delete [] T; 
    delete [] upperBoundary; 
    delete [] lowerBoundary; 
   
    // determine the program runtime in seconds 
    Real trun=(clock()-tstart)/CLOCKS_PER_SEC; 
    
    fprintf(dataFile,"nx = %d\n",nx); 
    fprintf(dataFile,"trun = %15.8f\n",trun); 
    fprintf(dataFile,"Tavg = %15.8f\n",sum); 
    fclose(dataFile); 

    // report successful run execution to console and return
    printf("Average temperature at t=%1.8f: %1.8f\n",t,sum); 
    printf("Runtime for nx=%d: %1.8f seconds\n",nx,trun); 

    return 0; 
} 

