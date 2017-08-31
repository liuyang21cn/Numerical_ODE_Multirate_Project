/*
 
 Yang Liu
 SMU Mathematcis
 
 2nd order SSP Partitioned Runge-Kutta method with m = 2
 
 */

#include "mat.h"
#include "MPRK2.h"
#include "rhs.h"

using namespace std;

// Assemble function is to put slow and fast parts into y
int MPRK2Stepper::Assemble(Mat &ynew, Mat &yslow, Mat &yfast, int const beg, int const end){
    
    int ntotal = yslow.Size()+yfast.Size();
    int nslow = yslow.Size();
    int nfast = yfast.Size();
    
    // beg and end are the begin point and end point of fast conponents with indices starts from 0
    if (beg == end){
        for (int i = 0; i<nslow; i++){
            ynew(i) = yslow(i);
        }
        ynew(nslow) = yfast(0);
    }
    else{
        
        for (int i = 0; i<beg; i++){ // 0-149
            ynew(i) = yslow(i);
        }
        for (int i = beg; i<=end; i++){
            ynew(i) = yfast(i-beg);
        }
        for (int i = end+1; i<ntotal; i++){
            ynew(i) = yslow(i-nfast);
        }
    }
    
    return 0;
}

// Split function is to pass information of y to slow and fast parts
int MPRK2Stepper::Split(Mat &y, Mat &yslow, Mat &yfast, int const beg, int const end){
    
    int ntotal = yslow.Size()+yfast.Size();
    int nslow = yslow.Size();
    int nfast = yfast.Size();
    
    if (beg == end){
        for (int i = 0; i<yslow.Size(); i++){
            yslow(i) = y(i);
        }
        yfast(0) = y(y.Size()-1);
    }
    else{
        
        for (int i = 0; i<beg; i++){
            yslow(i) = y(i);
        }
        for (int i = beg; i<=end; i++){
            yfast(i-beg) = y(i);
        }
        for (int i = end+1; i<ntotal; i++){
            yslow(i-nfast) = y(i);
        }
    }
    return 0;
    
}


Mat MPRK2Stepper::Evolve(Mat tspan, double h, Mat &y, int m, int beg, int end) {
    
    // check for legal inputs
    if (h <= 0.0) {
        std::cerr << "Multirate Base Method: Illegal h\n";
        return tspan;
    }
    if (tspan(1) <= tspan(0)) {
        std::cerr << "Multirate Base Method: Illegal tspan\n";
        return tspan;
    }
    
    // figure out how many time steps
    
    long int N = (tspan(1)-tspan(0))/h;
    if (tspan(1) > tspan(0)+N*h)  N++;
    
    // create ouput Mat
    Mat times(1, N+1);
    times(0) = tspan(0);
    
    // beg and end use C++ format indices, as 0 is the first index of array
    long int nfast = end-beg+1;
    long int nslow = y.Size()-nfast;

    // create local data
    Mat zs1(nslow), zs2(nslow), zf1(nfast), zf2(nfast);
    Mat ytemp(y.Size()), yn(y.Size()), y0(y);
    zs1 = 0.0;
    zs2 = 0.0;
    zf1 = 0.0;
    zf2 = 0.0;
    
    Mat fs1(nslow);
    Mat ff1(nfast);
    Mat fs2(nslow);
    Mat ff2(nfast);
   
    // iterate over time steps
    for (int i=0; i<N; i++) {
        
        // last step only: update h to stop directly at final time
        if (i == N-1)
            h = tspan(1)-times(i);
        
        // update y0 in every time step
        y0 = y;
        
        // iterate over all partitions
        for (int im = 0; im<m; im++) {
            
            // first part is little different
            if (im == 0) {
                
                // update zs1, zf1
                Split(y0, zs1, zf1, beg, end);
                
                // update f(zs1, zf1)
                if (frhs->Evaluate(times(i), y0, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }

                // update f(zs1) & f(zf1)
                Split(*f, fs1, ff1, beg, end);

                // update y at zs1, zf1
                y += 0.5*h/m*(*f);
                
                // update zs2, zf2
                zs2 = zs1 + h*fs1;
                zf2 = zf1 + h/m*ff1;
                
                // update f(zs2, zf2)
                Assemble(ytemp, zs2, zf2, beg, end);
                if (frhs->Evaluate(times(i)+(im+1.0)*h/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs2, ff2, beg, end);

                // update y at zs2, zf2
                y += 0.5*h/m*(*f);
                
            }
            else if (im > 0){
                
                // no need to update zs1
                
                // -- update zf1 ---
                Assemble(ytemp, zs1, zf1, beg, end);
                if (frhs->Evaluate(times(i)+(im-1.0)/m*h, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs1, ff1, beg, end);
                zf1 += 0.5*h/m*ff1;
                
                Assemble(ytemp, zs2, zf2, beg, end);
                if (frhs->Evaluate(times(i)+(im-0.5)/m*h, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs2, ff2, beg, end);
                zf1 += 0.5*h/m*ff2;
                
                // update f(zf1, zs1)
                Assemble(ytemp, zs1, zf1, beg, end);
                if (frhs->Evaluate(times(i)+im*h/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }

                // update y at zs1, zf1
                y += 0.5*h/m*(*f);
                
                // update zs2
                Assemble(ytemp, zs1, zf1, beg, end);
                if (frhs->Evaluate(times(i), ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs1, ff1, beg, end);
                zs2 = zs1 + h*fs1;
                
                // update zf2
                if (frhs->Evaluate(times(i)+h*im/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs1, ff1, beg, end);
                zf2 = zf1 + h/m*ff1;
                
                // update f(zs2, zf2)
                Assemble(ytemp, zs2, zf2, beg, end);
                if (frhs->Evaluate(times(i)+(im+1.0)*h/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                
                // update y at zs2, zf2
                y += 0.5*h/m*(*f);
                
            }
        }
        
        // update time
        times(i+1) = times(i) + h;
    }
    return times;
    
}


