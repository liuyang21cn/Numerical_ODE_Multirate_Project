/*
 
 Yang Liu
 SMU Mathematcis
 
 MPRK3 method
 */

#include "mat.h"
#include "MPRK3.h"
#include "rhs.h"

using namespace std;

// Assemble function is to put slow and fast parts into y
int MPRK3Stepper::Assemble(Mat &ynew, Mat &yslow, Mat &yfast, int const beg, int const end){
    
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
int MPRK3Stepper::Split(Mat &y, Mat &yslow, Mat &yfast, int const beg, int const end){
    
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

Mat MPRK3Stepper::Evolve(Mat tspan, double h, Mat &y, int m, const int beg, const int end) {
    
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
    Mat zs1(nslow), zs2(nslow), zs3(nslow), zs4(nslow);
    Mat zf1(nfast), zf2(nfast), zf3(nfast), zf4(nfast);
    Mat temps(nslow), tempf(nfast);
    
    Mat z1(nslow), z2(nfast), z3(nslow), z4(nfast), z5(nslow), z6(nfast), z7(nslow), z8(nfast);
    
    Mat ytemp(y.Size()), yn(y.Size()), y0(y), y0slow(nslow), y0fast(nfast);
    Mat yslow(nslow), yfast(nfast);
    zs1 = 0.0;
    zs2 = 0.0;
    zs3 = 0.0;
    zs4 = 0.0;
    zf1 = 0.0;
    zf2 = 0.0;
    zf3 = 0.0;
    zf4 = 0.0;
    
    Mat fs1(nslow), fs2(nslow), fs3(nslow), fs4(nslow);
    Mat ff1(nfast), ff2(nfast), ff3(nfast), ff4(nfast);

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
                Split(y0, y0slow, y0fast, beg, end);
                Split(y0, yslow, yfast, beg, end);
                
                // update f(zs1, zf1)
                if (frhs->Evaluate(times(i)+h*(im)/m, y0, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }

                // update f(zs1) & f(zf1)
                Split(*f, fs1, ff1, beg, end);
                
                // update y at zs1, zf1
                y += h*(1.0/6.0/m)*(*f);
                Split(y, yslow, yfast, beg, end);

                // update zs2, zf2
                zs2 = zs1 + h*(1.0/2.0)*fs1;
                zf2 = zf1 + h*(1.0/2.0/m)*ff1;
                
                // update f(zs2, zf2)
                Assemble(ytemp, zs2, zf2, beg, end);
                if (frhs->Evaluate(times(i)+h*(im+0.5)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff2, beg, end);
                
                if (frhs->Evaluate(times(i)+h*0.5, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs2, tempf, beg, end);
                
                // update yfast and yslow
                yfast += h*(1.0/3.0/m)*ff2;
                yslow += h*(1.0/3.0/m)*fs2;
                
                // update zs3
                zs3 = zs1 + h*(-1.0/6.0)*fs1 + h*(2.0/3.0)*fs2;
                
                // update zf3
                zf3 = zf1 + h*(-1.0/6.0/m)*ff1 + h*(2.0/3.0/m)*ff2;
                
                // update f(zs3, zf3)
                Assemble(ytemp, zs3, zf3, beg, end);
                if (frhs->Evaluate(times(i)+h*(im+0.5)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff3, beg, end);
                
                if (frhs->Evaluate(times(i)+h*0.5, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs3, tempf, beg, end);

                // update y at zs3, zf3
                yfast += h*(1.0/3.0/m)*ff3;
                yslow += h*(1.0/3.0/m)*fs3;
                
                // update zs4
                zs4 = zs1 + h*(1.0/3.0)*fs1 + h*(-1.0/3.0)*fs2 + h*1.0*fs3;

                // update zf4
                zf4 = zf1 + h*(1.0/3.0/m)*ff1 + h*(-1.0/3.0/m)*ff2 + h*(1.0/m)*ff3;
                
                // update f(zs4, zf4)
                Assemble(ytemp, zs4, zf4, beg, end);
                if (frhs->Evaluate(times(i)+h*(im+1.0)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff4, beg, end);
                
                if (frhs->Evaluate(times(i)+h, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs4, tempf, beg, end);
                
                // update y at zs3, zf3
                yfast += h*(1.0/6.0/m)*ff4;
                yslow += h*(1.0/6.0/m)*fs4;
                
                 Assemble(y, yslow, yfast, beg, end);
            }
            else if (im > 0){
                
                // -- update zs1 ---
                zs1 = y0slow;
                
                // -- update zf1 ---
                zf1 += h*(1.0/6.0/m)*ff1 + h*(1.0/3.0/m)*ff2 + h*(1.0/3.0/m)*ff3 + h*(1.0/6.0/m)*ff4;
                
                // update f(zs1, zf1)
                Assemble(ytemp, zs1, zf1, beg, end);
                if (frhs->Evaluate(times(i)+h*(im)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff1, beg, end);
                
                if (frhs->Evaluate(times(i), ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs1, tempf, beg, end);
                
                // update yfast and yslow
                yfast += h*(1.0/6.0/m)*ff1;
                yslow += h*(1.0/6.0/m)*fs1;
                
                // update zs2, zf2
                zs2 = zs1 + h*(1.0/2.0)*fs1;
                zf2 = zf1 + h*(1.0/2.0/m)*ff1;
                
                // update f(zs2, zf2)
                Assemble(ytemp, zs2, zf2, beg, end);
                if (frhs->Evaluate(times(i)+h*(im+0.5)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff2, beg, end);
                
                if (frhs->Evaluate(times(i)+h*0.5, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs2, tempf, beg, end);
                
                // update yfast and yslow
                yfast += h*(1.0/3.0/m)*ff2;
                yslow += h*(1.0/3.0/m)*fs2;
                
                // update zs3
                zs3 = zs1 + h*(-1.0/6.0)*fs1 + h*(2.0/3.0)*fs2;
                
                // update zf3
                zf3 = zf1 + h*(-1.0/6.0/m)*ff1 + h*(2.0/3.0/m)*ff2;
                
                // update f(zs3, zf3)
                Assemble(ytemp, zs3, zf3, beg, end);
                if (frhs->Evaluate(times(i)+h*(im+0.5)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff3, beg, end);
                
                if (frhs->Evaluate(times(i)+h*0.5, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs3, tempf, beg, end);
                
                // update y at zs3, zf3
                yfast += h*(1.0/3.0/m)*ff3;
                yslow += h*(1.0/3.0/m)*fs3;
                
                // update zs4
                zs4 = zs1 + h*(1.0/3.0)*fs1 + h*(-1.0/3.0)*fs2 + h*1.0*fs3;
                
                // update zf4
                zf4 = zf1 + h*(1.0/3.0/m)*ff1 + h*(-1.0/3.0/m)*ff2 + h*(1.0/m)*ff3;
                
                // update f(zs4, zf4)
                Assemble(ytemp, zs4, zf4, beg, end);
                if (frhs->Evaluate(times(i)+h*(im+1.0)/m, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, temps, ff4, beg, end);
                
                if (frhs->Evaluate(times(i)+h, ytemp, *f) != 0) {
                    std::cerr << "Multirate Partitioned RK2: Error in 1st ODE RHS function\n";
                    return times;
                }
                Split(*f, fs4, tempf, beg, end);
                
                // update y at zs3, zf3
                yfast += h*(1.0/6.0/m)*ff4;
                yslow += h*(1.0/6.0/m)*fs4;
                
                Assemble(y, yslow, yfast, beg, end);
            }
        }
        
        // update time
        times(i+1) = times(i) + h;
    }
    return times;
    
}


