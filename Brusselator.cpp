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
#include "erk4.h"


using namespace std;

//----- Brusselator equations -----//
class RHS: public RHSFunction {
    
public:
    double a, b, epsilon;
    double u, v, w;
    int Evaluate(double x, Mat &y, Mat &f) {
        
        u = y(0);
        v = y(1);
        w = y(2);
        
        f(0) = a - (w+1.0)*u + v*u*u;
        f(1) = w*u - v*u*u;
        f(2) = (b-w)/epsilon - w*u;
        
        return 0;
    }
};
// ---------------------------------------------//


// main routine
int main() {
    
    // time steps to try
    Mat h(" 0.01, 0.005, 0.0025, 0.00125");
    
    // step for fast component
    Mat m("1, 2, 4, 16, 64, 256, 512, 1024");
    
    //----- brusselator equations -----//
    // initial condition
    Mat y0(3);
    y0(0) = 3.9;
    y0(1) = 1.1;
    y0(2) = 2.8;
    
    // set up step
    double t0 = 0.0;
    double Tf = 10.0;
    double tcur = t0;
    double dtout = 0.1;
    
    // true solution
    Mat yt(3);
    
    // create ODE RHS function objects
    RHS rhs;
    rhs.a = 1.2;
    rhs.b = 2.5;
    rhs.epsilon = 0.01;
    
    // create multirate methdod objects
    MultirateBaseMethod MBM(rhs, y0);
    MPRK2Stepper        MPRK2(rhs, y0);
    MPRK3Stepper        MPRK3(rhs, y0);
    ERK4Stepper         ERK4(rhs);

    // temporary variables
    Mat y1(3), y2(3), y3(3), yerr1(3), yerr2(3), yerr3(3), tspan(1,2);
    double err, maxerr;
    Mat errs(h.Size());
    //-----------------------------------------------//

    // true solution //
    // loop over time step sizes
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            yt = y0;
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
                Mat tvals = ERK4.Evolve(tspan, h(ih), yt);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
            }
            
        }
    
    //---------------------------------//
    
   cout << "\nA system with a slow variable coupled into the fast equations:";
    
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
                Mat tvals = MBM.Evolve(tspan, h(ih), y1, m(im), 2, 2);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                

            }
            // compute the error at tcur, output to screen and accumulate maximum
            yerr1 = y1 - yt;
            err = yerr1.MaxNorm();
            maxerr = (maxerr > err) ? maxerr : err;
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
                Mat tvals = MPRK2.Evolve(tspan, h(ih), y2, m(im), 2, 2);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                

            }
            // compute the error at tcur, output to screen and accumulate maximum
            yerr2 = y2 - yt;
            err = yerr2.MaxNorm();
            maxerr = (maxerr > err) ? maxerr : err;
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
                Mat tvals = MPRK3.Evolve(tspan, h(ih), y3, m(im), 2, 2);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                

            }
            // compute the error at tcur, output to screen and accumulate maximum
            yerr3 = y3 - yt;
            err = yerr3.MaxNorm();
            maxerr = (maxerr > err) ? maxerr : err;
            cout << "   h = " << h(ih) << "\t  error = " << maxerr;
            errs(ih) = maxerr;
            if (ih > 0)
                cout << "\t  conv rate = " << (log(errs(ih))-log(errs(ih-1)))/(log(h(ih))-log(h(ih-1)));
            cout << endl;
        }
    }
    
    return 0;
}
