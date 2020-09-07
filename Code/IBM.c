/**
 * SIMULATE FROM THE IBM
 *
 * author: Alexander P. Browning (ap.browning@qut.edu.au)
 *         School of Mathematical Sciences 
 *         Queensland University of Technology
 *
 * Date: 18/10/2019
 */

#include "mex.h"
#include <math.h>
#include "matrix.h"

// Universal constants
#define pi          3.14159265358979323846
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define max(a,b)    (((a) < (b)) ? (b) : (a))

/* SETTINGS */

// Max agents
const int Nmax = 5000;

// Density profile settings
const int nBins = 80;
double BinWidth = 23.75;

// Pair correlation settings
double PC_dr        = 5;
const int PC_nBins  = 20;

// Output interval
double dtout        = 0.5;
const int numTout   = 73;

/* END SETTINGS */

/**
 * RANDOM NUMBER GENERATOR
 */
int rand2();
int rseed = 0;
inline void srand2(int x) {
    rseed = x;
}
#define RAND_MAX2 ((1U << 31) - 1)
inline int rand2() {
    return rseed = (rseed * 1103515245 + 12345) & RAND_MAX2;
}



/**
 * UNIFORM NUMBER GENERATOR
 */
double sampleU() {
    double random = ((double)rand2()) / ((double)RAND_MAX2);
    return random;
}


/**
 * VM SAMPLE
 */
double VM(double x, double mu, double kappa) {

    double out = exp(kappa * cos(x - mu));
    return out;

}
double sampleVM(double mu, double kappa) {

    // Sample from non-normalised density with rejection sampling
    double fmax = VM(mu,mu,kappa);

    double x = sampleU()*2*pi - pi;
    double y = VM(x,mu,kappa);
    double u = sampleU()*fmax;

    while (u > y) {
        x = sampleU()*2*pi - pi;
        y = VM(x,mu,kappa);
        u = sampleU()*fmax;
    }

    return x;

}


/**
 * PERIODIC SQUARE DISTANCE BETWEEN AGENTS
 */
double distance2(double x1, double x2, double y1, double y2, double L, double H) {

    double dx = fabs(x1 - x2);
    double dy = fabs(y1 - y2);

    dx = min(dx,L - dx);
    dy = min(dy,H - dy);

    double d = pow(dx,2) + pow(dy,2);
    return d;

}


/**
 * 1D PERIODIC DISPLACEMENT (x1 -> x2)
 */
double disp(double x1, double x2, double L) {

    double s = 0;

    double dx = x2 - x1;
    double adx = fabs(dx);
    double dxL = L - adx;

    if (adx < dxL) {
        s = dx;
    } else {
        if (x1 < L / 2) {
            s = -dxL;
        } else {
            s = dxL;
        }
    }

    return s;

}


/**
 * GAUSSIAN KERNEL
 */
double kernel(double r2, double sigma2, double gamma) {

    double y = 0;

    if (r2 < 9 * sigma2) {

        y = gamma * exp( -r2 / (2 * sigma2) );

    }

    return y;

}


/*
 * RANDOMLY CHOOSE
 */
int choose_agent(double* M, double M_tot) {

    int i = 0;
    double Mc = max(0,M[i]);
    double alpha = sampleU() * M_tot;
    while (alpha > Mc) {
        i += 1;
        Mc += max(0,M[i]);
    }
    return i;

}


/* 
 * PERIODIC MODULUS
 */
double mod(double x, double L) {

    if (x <= 0) {
        x += L;
    } else if (x >= L) {
        x -= L;
    }

    return x;

}


/*
 * CALCULATE PAIR CORRELATION
 *  AVERAGE OF PERODIC PC BETWEEN Y < 400 AND Y > 1500
 */
int PairCorrelation(int N, double * X, double * Y, double L, double * PC) {
    
    // Initialise
    double PC1[PC_nBins];
    double PC2[PC_nBins];
    for (int i = 0; i < PC_nBins; i++) {
        PC1[i] = 0;
        PC2[i] = 0;
    }
    
    double x0, y0, x1, y1, dx, dy, dist;

    int N1 = 0;      // N < 400
    int N2 = 0;      // N > 1500
    
    
    // Loop through agents 1 (x0,y0)
    for (int i = 0; i < N; i++) {
        x0 = X[i];
        y0 = Y[i];
        
        // y0 < 400 (and only look at second agent with y1 < 400)
        if (y0 < 400) {
            N1 += 1;
            
            // Loop through agents 2 (x1,y1)
            for (int j = 0; j < N; j++) { if (i != j) {
                double x1 = X[j];
                double y1 = Y[j];
                
                if (y1 < 400) {
                    
                    dx = fabs(x1 - x0);
                    dx = min(dx,L - dx);
                    
                    dy = fabs(y1 - y0);
                    dy = min(dy,400 - dy);
                    
                    dist = sqrt(dx * dx + dy * dy);
                    
                    // Bin distance
                    for (int k = 0; k < PC_nBins; k++) {
                        if (dist >= k * PC_dr && dist < (k+1) * PC_dr) {
                            PC1[k] += 1;
                            break;
                        }
                    }
                    
                } // end (y1 < 400)
                
            }} // end agent loop 2
                        
        } // end agent loop 1 (y0 < 400)
        
               
        // y0 > 1500 (and only look at second agent with y1 > 1500)
        if (y0 > 1500) {
            N2 += 1;
            
            // Loop through agents 2 (x1,y1)
            for (int j = 0; j < N; j++) { if (i != j) {
                double x1 = X[j];
                double y1 = Y[j];
                
                if (y1 > 1500) {
                    
                    dx = fabs(x1 - x0);
                    dx = min(dx,L - dx);
                    
                    dy = fabs(y1 - y0);
                    dy = min(dy,400 - dy);
                    
                    dist = sqrt(dx * dx + dy * dy);
                    
                    // Bin distance
                    for (int k = 0; k < PC_nBins; k++) {
                        if (dist >= k * PC_dr && dist < (k+1) * PC_dr) {
                            PC2[k] += 1;
                            break;
                        }
                    }
                    
                } // end (y1 > 1500)
                
            }} // end agent loop 2
                        
        } // end agent loop 1 (y0 > 1500)
           
    }    
    
    
    // Average and scale to get PC
    for (int i = 0; i < PC_nBins; i++) {

        PC1[i] /= N1 * N1 / (L * 400) * pi * (pow((i+1) * PC_dr,2) - pow(i * PC_dr,2));
        PC2[i] /= N2 * N2 / (L * 400) * pi * (pow((i+1) * PC_dr,2) - pow(i * PC_dr,2));
        
        PC[i] = 0.5 * (PC1[i] + PC2[i]);
        
    }
    
    
    return 1;
    
}


/*
 * CALCULATE DENSITY PROFILE
 */
int DensityProfile(int N, double * Y, double * D) {
    
    double BinStart;
    double yloc;
    
    // Start D
    for (int i = 0; i < nBins; i++) {
        D[i] = 0;
    }
    
    // Loop through agents
    for (int agent = 0; agent < N; agent++) {
        yloc = Y[agent];
        
        // Determine appropriate bin
        for (int i = 0; i < nBins; i++) {
            BinStart = i * BinWidth;
            if (yloc > BinStart && yloc <= (BinStart + BinWidth)) {
                D[i] += 1;
            }
        }
               
    }
    
    return 1;
    
}
    


/**
 * SIMULATE BINNY MODEL
 */
int Binny(double* params, double* domain, double* IC, int N0, double T, double Nmaxsoft, double * X, double * Y, int * N18, double * X18, double * Y18, int * NTout)
{
    
    // GET PARAMETERS
    double m        = *(params);
    double p        = *(params + 1);
    double gm       = *(params + 2);
    double gp       = *(params + 3);
    double gb       = *(params + 4);
    double s2       = pow(*(params + 5),2); // default: 24
    double mu_s     = *(params + 6); // default: 24
    
    double L        = *(domain);
    double H        = *(domain + 1);
    
    // END GET PARAMETERS
    
    // INITIALISE
    int     N       = N0;
    double  t       = 0;
    
    NTout[0]        = N0;
    int next_tout   = 1;
    
    double M[Nmax];
    double P[Nmax];
    double Bx[Nmax];
    double By[Nmax];
    
    // t = 18h OUTPUT FLAG
    int output18 = 0;
    
    // Loop through each agent
    for (int i = 0; i < N0; i++) {
        
        double x1 = *(IC + i);
        double y1 = *(IC + N0 + i);
        
        double MrS = m;
        double PrS = p;
        double BxS = 0;
        double ByS = 0;
        
        // Loop through other agents
        for (int j = 0; j < N; j++) { if (i != j) {
            
            double x2 = *(IC + j);
            double y2 = *(IC + N0 + j);
            
            double r2 = distance2(x1,x2,y1,y2,L,H);
            
            double b_ = kernel(r2,s2,gb) / s2;
            MrS      -= kernel(r2,s2,gm);
            PrS      -= kernel(r2,s2,gp);
            
            if (b_ != 0) {
                // Bx is b_ * (disp from them to us)
                BxS += b_ * disp(x2,x1,L);
                ByS += b_ * disp(y2,y1,H);
            }
            
        }}
                
        X[i] = x1;
        Y[i] = y1;
        
        M[i] = MrS;
        P[i] = PrS;
        Bx[i] = BxS;
        By[i] = ByS;
        
    }

    // END INITIALISE

        
    
    
    // LOOP THROUGH TIME
    while (t < T && N < Nmax && N < Nmaxsoft) {
        
        // CALCULATE TOTAL EVENT RATES
        double M_tot = 0;
        double P_tot = 0;
        for (int i = 0; i < N; i++) {    
            M_tot += max(0,M[i]);
            P_tot += max(0,P[i]);
        }
       
        // SAMPLE TIMESTEP
        double tau = -log(sampleU()) / (M_tot + P_tot);
        t += tau;
        
        // STOP IF NEXT EVENT OCCURS AFTER t = T
        if (t > T) {
            break;
        }
        
        // DECIDE EVENT
        double alpha = sampleU() * (M_tot + P_tot);
        
        // MOVEMENT
        if (alpha < M_tot) {
                     
            // CHOOSE AGENT
            int i = choose_agent(M,M_tot);
            
            // LOCATION
            double xc = X[i];
            double yc = Y[i];
            
            // INCREASE RATES OF SURROUNDING AGENTS
            for (int j = 0; j < N; j++) { if (i != j) {
                
                double x2 = X[j];
                double y2 = Y[j];
                
                double r2 = distance2(xc,x2,yc,y2,L,H);
                double m_ = kernel(r2,s2,gm);
                double p_ = kernel(r2,s2,gp);
                double b_ = kernel(r2,s2,gb) / s2;
                
                if (m_ != 0) {
                    M[j] += m_;
                }
                if (p_ != 0) {
                    P[j] += p_;
                }
                if (b_ != 0) {
                    // Bx is b_ * (disp from them to us)
                    Bx[j] -= b_ * disp(xc,x2,L);
                    By[j] -= b_ * disp(yc,y2,H);
                }
                
            }}
            // END INCREASE RATES
            
            
            // MOVE SOMEWHERE, INCLUDE BIAS
            double md       = mu_s;
            double Bx_i     = Bx[i];
            double By_i     = By[i];
            double vm_mu    = atan2(By_i,Bx_i);
            double vm_kappa = sqrt(pow(Bx_i,2) + pow(By_i,2));
            
            double theta    = sampleVM(vm_mu,vm_kappa);
            
            double xp       = mod(xc + md * cos(theta),L);
            double yp       = mod(yc + md * sin(theta),H);
            
            // UPDATE RATES OF AGENT
            double MrS      = m;
            double PrS      = p;
            double BxS      = 0;
            double ByS      = 0;
            
            for (int j = 0; j < N; j++) { if(i != j) {

                double x2 = X[j];
                double y2 = Y[j];

                double r2 = distance2(xp,x2,yp,y2,L,H);
                double m_ = kernel(r2,s2,gm);
                double p_ = kernel(r2,s2,gp);
                double b_ = kernel(r2,s2,gb) / s2;

                if (m_ != 0) {
                    M[j] -= m_;
                    MrS  -= m_;
                }
                if (p_ != 0) {
                    P[j] -= p_;
                    PrS  -= p_;
                }
                if (b_ != 0) {
                    // Bx is b_ * (disp from them to us)
                    double sx = disp(x2,xp,L);
                    double sy = disp(y2,yp,H);
                    
                    // This agents bias
                    BxS      += b_ * sx;
                    ByS      += b_ * sy;
                    
                    // Other agents bias (displacement -ve)
                    Bx[j]    -= b_ * sx;
                    By[j]    -= b_ * sy;
                }
               
            }}
            // END UPDATE RATES
            
            
            // "MOVE" AGENT
            X[i]    = xp;
            Y[i]    = yp;
            M[i]    = MrS;
            P[i]    = PrS;
            Bx[i]   = BxS;
            By[i]   = ByS;
            
            
       // PROLIFERATION
       } else {
            
            // CHOOSE AGENT
            int i = choose_agent(P,P_tot);
            
            // LOCATION
            double xc = X[i];
            double yc = Y[i];
            
            // NEW LOCATION (USE BIAS)
            double md       = mu_s;
            double Bx_i     = Bx[i];
            double By_i     = By[i];
            double vm_mu    = atan2(By_i,Bx_i);
            double vm_kappa = sqrt(pow(Bx_i,2) + pow(By_i,2));
            
            double theta    = sampleVM(vm_mu,vm_kappa);
            
            double xp       = mod(xc + md * cos(theta),L);
            double yp       = mod(yc + md * sin(theta),H);
            
            // OLD (BIVARIATE NORMAL)
            //double u1 = sampleU();
            //double u2 = sampleU();
            //
            //double xp = mod(xc + sigma_d * sqrt(-2 * log(u1)) * cos(2*pi*u2),L);
            //double yp = mod(yc + sigma_d * sqrt(-2 * log(u1)) * sin(2*pi*u2),H);           
            
            
            // UPDATE RATES
            double MrS = m;
            double PrS = p;
            double BxS = 0;
            double ByS = 0;
            
            for (int j = 0; j < N; j++) {
             
                double x2 = X[j];
                double y2 = Y[j];

                double r2 = distance2(xp,x2,yp,y2,L,H);
                double m_ = kernel(r2,s2,gm);
                double p_ = kernel(r2,s2,gp);
                double b_ = kernel(r2,s2,gb) / s2;

                if (m_ != 0) {
                    M[j] -= m_;
                    MrS  -= m_;
                }
                if (p_ != 0) {
                    P[j] -= p_;
                    PrS  -= p_;
                }
                if (b_ != 0) {

                    // Bx is b_ * (disp from them to us)
                    double sx = disp(x2,xp,L);
                    double sy = disp(y2,yp,H);
                    
                    // This agents bias
                    BxS      += b_ * sx;
                    ByS      += b_ * sy;
                    
                    // Other agents bias (displacement -ve)
                    Bx[j]    -= b_ * sx;
                    By[j]    -= b_ * sy;
                    
                }
                
            }
            // END UPDATE RATES
            
            // "CREATE" NEW AGENT
             X[N]    = xp;
             Y[N]    = yp;
             M[N]    = MrS;
             P[N]    = PrS;
             Bx[N]   = BxS;
             By[N]   = ByS;
             N += 1;
            
        } // END PROLIFERATION
        
        // OUTPUT MID-TIME DATA
        if (t > 18 && output18 == 0) {
            
            N18[0] = N;
            for (int i = 0; i < N; i++) {
                
                X18[i] = X[i];
                Y18[i] = Y[i];
                        
            } 
            output18 = 1;
            
        }
        
        // OUTPUT TIME
        if (t > next_tout * dtout) {
            
            NTout[next_tout] = N;
            next_tout += 1;
            
        }

    }

    NTout[numTout-1] = N;
    
    return N;
    
}
// END BINNY FUNCTION


/* gateway function
 *
 * This function gets input from MATLAB and returns output to MATLAB.
 *
 * Inputs:
 *  params  : [m,p,gm,gp,gb,sig,mu_s]           : (1 x 7)    vector
 *  domain  : [L,H]                             : (1 x 2)    vector
 *  IC      : [X1,Y1; ...; XN,YN]               : (N0 x 2) matrix
 *  T       : (36)                              : Solve for 0 < t < T
 *  Nmax    : (5000)                            : Maximum population before cutout
 *  seed    : (int)                             : rng seed
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("mex:nrhs","Requires six inputs: params,domain,IC,T,Nmax,seed");
    }

    if (nlhs != 1 && nlhs != 2 && nlhs != 6) {
        mexErrMsgIdAndTxt("mex:nlhs","Either one, two or six outputs required.");
    }

    /* Assumes user enters valid inputs*/
    double * inparams   = (double*) mxGetDoubles(prhs[0]);
    mwSize   ninparams  = (int)     mxGetN(prhs[0]);
    double * domain     = (double*) mxGetDoubles(prhs[1]);
    double * IC         = (double*) mxGetDoubles(prhs[2]);
    mwSize   N0         = (int)     mxGetM(prhs[2]);
    double   T          = (double)  mxGetScalar(prhs[3]);
    int      Nmaxsoft   = (int)     mxGetScalar(prhs[4]);
    int      seed       = (int)     mxGetScalar(prhs[5]);

    // Process input parameters
    if (ninparams != 7) {
        mexErrMsgIdAndTxt("mex:nrhs","params must have a length of seven.");
    }
    double params[7];
    for (int i = 0; i < ninparams; i++) {
        params[i] = inparams[i];
    }
    
    // Random seed
    srand2(seed);
    
    // Simulate model
    double X[Nmax];
    double Y[Nmax];
    int    N18[1];
    double X18[Nmax];
    double Y18[Nmax];
    int NTout[numTout];
    int NT          = Binny(params,domain,IC,N0,T,Nmaxsoft,X,Y,N18,X18,Y18,NTout);
        
    // One outout: Nt
    if (nlhs == 1) {
        
        // Create output matrix
        plhs[0] =  mxCreateDoubleMatrix(numTout,1,mxREAL);
        double * NToutML = mxGetDoubles(plhs[0]);
    
        // Loop through time points
        for (int i = 0; i < numTout; i++) {
     
            NToutML[i] = NTout[i];
        
        }
        
    }
    
    // Two output: Locations
    if (nlhs == 2) {
        
        // Create output matrix
        plhs[0]         = mxCreateDoubleMatrix(NT,2,mxREAL);
        plhs[1]         = mxCreateDoubleMatrix(N18[0],2,mxREAL);
        double * XY     = mxGetDoubles(plhs[0]);
        double * XY18   = mxGetDoubles(plhs[1]);

        // Loop through agents to output
        for (int i = 0; i < NT; i++) {

            // X positions
            XY[i]    = X[i];
            XY[NT+i] = Y[i];

        }
        
        // Loop through agents to output
        for (int i = 0; i < N18[0]; i++) {

            // X positions
            XY18[i]    = X18[i];
            XY18[N18[0]+i] = Y18[i];

        }

    }
    
    // Three outputs: [N,PC,D]
    if (nlhs == 6) {
        
        // Initialise outputs
        plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
        plhs[1] = mxCreateDoubleMatrix(PC_nBins,1,mxREAL);
        plhs[2] = mxCreateDoubleMatrix(nBins,1,mxREAL);
        
        plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
        plhs[4] = mxCreateDoubleMatrix(PC_nBins,1,mxREAL);
        plhs[5] = mxCreateDoubleMatrix(nBins,1,mxREAL);
                
        double * Nout   = mxGetDoubles(plhs[0]);
        double * PCout  = mxGetDoubles(plhs[1]);
        double * Dout   = mxGetDoubles(plhs[2]);
        
        double * Nout18   = mxGetDoubles(plhs[3]);
        double * PCout18  = mxGetDoubles(plhs[4]);
        double * Dout18   = mxGetDoubles(plhs[5]);
        
        // Fill outputs
        Nout[0]         = NT;
        DensityProfile(NT,Y,Dout);
        PairCorrelation(NT, X, Y, *(domain), PCout);
        
        Nout18[0]       = N18[0];
        DensityProfile(N18[0],Y18,Dout18);
        PairCorrelation(N18[0],X18,Y18, *(domain), PCout18);
        
    }
    
    
}
