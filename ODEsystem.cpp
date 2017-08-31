/*
 
 Yang Liu
 SMU Mathematics
 
 test file for multirate project
 
 */

#include <iostream>
#include <math.h>
#include "mat.h"
#include "MBM.h"
#include "MPRK2.h"
#include "MPRK3.h"

using namespace std;

//----- A one-way system with the fast variables coupled into the slow equation -----//
class RHS: public RHSFunction {
    
public:
    double x, y, z;
    int Evaluate(double t, Mat &yy, Mat &f) {
        x = yy(0);
        y = yy(1);
        z = yy(2);
        
        f(0) = -50*y;
        f(1) = 50*x;
        f(2) = -z+x+y;
        
        return 0;
    }
};

// Exact solution of Prothero-Robinson Problem
Mat ytrue(const double t, const double omega) {
    Mat yt(3);
    yt(0) = cos(50.0*t);
    yt(1) = sin(50.0*t);
    yt(2) = 5051.0/2501.0*exp(-t)-49.0/2501.0*cos(50.0*t)+51.0/2501.0*sin(50.0*t);
    return yt;
};
// ---------------------------------------------//


// main routine
int main() {
    
    // time steps to try
    Mat h(" 0.01, 0.005, 0.0025, 0.00125");
    
    // step for fast component
    Mat m("1, 2, 4, 16, 64, 256, 512, 1024");
    
    //----- A one-way system with the fast variables coupled into the slow equation -----//
    // initial condition
    Mat y0(3);
    y0(0) = 1;
    y0(1) = 0;
    y0(2) = 2;
    
    // set up step
    double t0 = 0.0;
    double Tf = 1.0;
    double tcur = t0;
    double dtout = 0.1;
    
    // true solution
    Mat yt(3);
    
    // create ODE RHS function objects
    RHS rhs;
    
    // create multirate bset methdod objects
    MultirateBaseMethod MBM(rhs, y0);
    MPRK2Stepper        MPRK2(rhs, y0);
    MPRK3Stepper        MPRK3(rhs, y0);
    
    
    // temporary variables
    Mat y1(3), y2(3), y3(3), y4(3), yerr1(3), yerr2(3), yerr3(3), yerr4(3), tspan(1,2);
    double err, maxerr;
    Mat errs(h.Size());
    //-----------------------------------------------//

    
    cout << "\nA one-way system with the fast variables coupled into the slow equation:";
    
    //----- Multirate Base Method -----//
    cout << "\n 1.Multirate Base Method :\n";

    // loop over time step sizes
    for (int im = 0; im<m.Size(); im++) {
        cout<< "\n m = " << m(im) << endl;
        
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            y1 = y0;
            tcur = t0;
            
            // reset maxerr
            maxerr = 0.0;
            
            // loop over output step sizes: call solver and output error
            while (tcur < 0.99999*Tf) {
                
                // set the time interval for this solve
                tspan(0) = tcur;
                tspan(1) = tcur + dtout;
                
                if (tspan(1) > Tf)  tspan(1) = Tf;
                
                // call the solver, update current time
                Mat tvals = MBM.Evolve(tspan, h(ih), y1, m(im), 0, 1);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
                // compute the error at tcur, output to screen and accumulate maximum
                yt = ytrue(tcur, 5.0);
                yerr1 = y1 - yt;
                err = yerr1.MaxNorm();
                maxerr = (maxerr > err) ? maxerr : err;
            }
            
            cout << "   h = " << h(ih) << "\t  error = " << maxerr;
            errs(ih) = maxerr;
            if (ih > 0)
                cout << "\t  conv rate = " << (log(errs(ih))-log(errs(ih-1)))/(log(h(ih))-log(h(ih-1)));
            cout << endl;
        }
    }
    //---------------------------------//
    
   //----- Constantinescu and Sandu RK2a Method ----------//
    cout << "\n 2.Constantinescu and Sandu RK2a Method :\n";

    // loop over time step sizes
    for (int im = 0; im<m.Size(); im++) {
        cout<< "\n m = " << m(im) << endl;
        
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            y2 = y0;
            tcur = t0;
            
            // reset maxerr
            maxerr = 0.0;
            
            // loop over output step sizes: call solver and output error
            while (tcur < 0.99999*Tf) {
                
                // set the time interval for this solve
                tspan(0) = tcur;
                tspan(1) = tcur + dtout;
                
                if (tspan(1) > Tf)  tspan(1) = Tf;
                
                // call the solver, update current time
                Mat tvals = MPRK2.Evolve(tspan, h(ih), y2, m(im), 0, 1);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
                // compute the error at tcur, output to screen and accumulate maximum
                yt = ytrue(tcur, 5.0);
                yerr2 = y2 - yt;
                err = yerr2.MaxNorm();
                maxerr = (maxerr > err) ? maxerr : err;
            }
            
            cout << "   h = " << h(ih) << "\t  error = " << maxerr;
            errs(ih) = maxerr;
            if (ih > 0)
                cout << "\t  conv rate = " << (log(errs(ih))-log(errs(ih-1)))/(log(h(ih))-log(h(ih-1)));
            cout << endl;
        }
    }
    // -----------------------------------------//
    
    
    //----- Constantinescu and Sandu RK43 Method ----------//
    cout << "\n 3.Constantinescu and Sandu RK43 :\n";
    
    // loop over time step sizes
    for (int im = 0; im<m.Size(); im++) {
        cout<< "\n m = " << m(im) << endl;
        
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            y3 = y0;
            tcur = t0;
            
            // reset maxerr
            maxerr = 0.0;
            
            // loop over output step sizes: call solver and output error
            while (tcur < 0.99999*Tf) {
                
                // set the time interval for this solve
                tspan(0) = tcur;
                tspan(1) = tcur + dtout;
                
                if (tspan(1) > Tf)  tspan(1) = Tf;
                
                // call the solver, update current time
                Mat tvals = MPRK3.Evolve(tspan, h(ih), y3, m(im), 0, 1);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
                // compute the error at tcur, output to screen and accumulate maximum
                yt = ytrue(tcur, 5.0);
                yerr3 = y3 - yt;
                err = yerr3.MaxNorm();
                maxerr = (maxerr > err) ? maxerr : err;
            }
            
            cout << "   h = " << h(ih) << "\t  error = " << maxerr;
            errs(ih) = maxerr;
            if (ih > 0)
                cout << "\t  conv rate = " << (log(errs(ih))-log(errs(ih-1)))/(log(h(ih))-log(h(ih-1)));
            cout << endl;
        }
    }
    
    return 0;
}
