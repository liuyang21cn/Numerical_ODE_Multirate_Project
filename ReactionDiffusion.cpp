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

//----- Advection Diffusion Problem -----//
class RHS_RD: public RHSFunction {
    
public:
    
    int Evaluate(double t, Mat &y, Mat &f) {
        
        double epsilon = 0.01;
        double gamma = 100;
        double lambda = 0.5*sqrt(2*gamma/epsilon);
        double pi = 3.14159265359;
        
        int Ny = y.Size();
        double dx = 3.5/(Ny-1);
  
        
        f = 0.0;
        for (int i = 1; i<Ny-1; i++){
            f(i) = epsilon/dx/dx*(y(i+1)-2*y(i)+y(i-1)) + gamma*y(i)*y(i)*(1-y(i));
        }
        
        f(0) = 0;
        f(Ny-1) = 0;
        
        return 0;
    }
};

// ---------------------------------------------//

// main routine
int main() {
    
    // time steps to try
    Mat h(" 0.01, 0.005, 0.001");
    
    // step for fast component
    //    Mat m("100, 200, 400");
    Mat m("1, 2, 4, 16, 64, 256");
    
    //----- Advection Diffusion Problem -----//
    // initial condition
    int Ny = 401;
    double L = 3.5;
    
    double epsilon = 0.01;
    double gamma = 100;
    double lambda = 0.5*sqrt(2*gamma/epsilon);
    double pi = 3.14159265359;
    
    double dx = 3.5/(Ny-1);
    double x = 0.0;
    
    Mat y0(Ny);
    for (int i = 0; i<Ny; i++){
        y0(i) = 1.0/(1.0+exp(lambda*(x-1.0)));
        x += dx;
    }
    
    // set up step
    double t0 = 0.0;
    double Tf = 0.5;
    double tcur = t0;
    double dtout = 0.1;
    
    // create ODE RHS function objects
    RHS_RD rhs;
    
    // create multirate methdod objects
    MultirateBaseMethod MBM(rhs, y0);
    MPRK2Stepper        MPRK2(rhs, y0);
    MPRK3Stepper        MPRK3(rhs, y0);
    
    // temporary variables
    Mat y1(Ny), y2(Ny), y3(Ny), y4(Ny), yerr1(Ny), yerr2(Ny), yerr3(Ny), yerr4(Ny), tspan(1,2);
    double err1, maxerr1;
    Mat errs1(h.Size()), errs2(h.Size());
    
    //-----------------------------------------//
    
    
    ///////// Advection-Diffusion equation /////////
    cout << "\nAdvection-Diffusion equation problem:\n";
    
    //----- Multirate Base Method -----//
    cout << "\n 1. Multirate Base Method:\n";
    
    // loop over time step sizes
    for (int im = 0; im<m.Size(); im++) {
        cout<< "\n m = " << m(im) << endl;
        
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            y1 = y0;
            tcur = t0;
            
            // reset maxerr
            maxerr1 = 0.0;
            
            // loop over output step sizes: call solver and output error
            while (tcur < 0.99999*Tf) {
            
                // set the time interval for this solve
                tspan(0) = tcur;
                tspan(1) = tcur + dtout;
                if (tspan(1) > Tf)  tspan(1) = Tf;
                
                // call the solver, update current time
                Mat tvals = MBM.Evolve(tspan, h(ih), y1, m(im), 10, 380);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals

            }
            //errs1(ih) = y1(149) -
            cout << "   h = " << h(ih);
            cout << "y1(150) = " << y1(149) << ", y1(155) = " << y1(154) << endl;
        }
    }
    //--------------------------------------------//
    
    //----- Constantinescu and Sandu RK2a Method -----//
    cout << "\n 2. Constantinescu and Sandu RK2a Method:\n";
    
    // loop over time step sizes
    for (int im = 0; im<m.Size(); im++) {
        cout<< "\n m = " << m(im) << endl;
        
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            y2 = y0;
            tcur = t0;
            
            // reset maxerr
            maxerr1 = 0.0;
            
            // loop over output step sizes: call solver and output error
            while (tcur < 0.99999*Tf) {
                
                // set the time interval for this solve
                tspan(0) = tcur;
                tspan(1) = tcur + dtout;
                if (tspan(1) > Tf)  tspan(1) = Tf;
                
                // call the solver, update current time
                Mat tvals = MPRK2.Evolve(tspan, h(ih), y2, m(im), 10, 380);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
            }
            cout << "   h = " << h(ih);

            cout << "y2(150) = " << y2(149) << ", y2(155) = " << y2(154) << endl;
        }
    }
    //----------------------------------------------------//
    
    //----- Constantinescu and Sandu RK43 Method -----//
    cout << "\n 3. Constantinescu and Sandu RK43 Method:\n";
    
    // loop over time step sizes
    for (int im = 0; im<m.Size(); im++) {
        cout<< "\n m = " << m(im) << endl;
        
        for (int ih=0; ih<h.Size(); ih++) {
            
            // set the initial condition, initial time
            y3 = y0;
            tcur = t0;
            
            // reset maxerr
            maxerr1 = 0.0;
            
            // loop over output step sizes: call solver and output error
            while (tcur < 0.99999*Tf) {
                
                // set the time interval for this solve
                tspan(0) = tcur;
                tspan(1) = tcur + dtout;
                if (tspan(1) > Tf)  tspan(1) = Tf;
                
                // call the solver, update current time
                Mat tvals = MPRK3.Evolve(tspan, h(ih), y3, m(im), 10, 380);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
            }
            cout << "   h = " << h(ih);

            cout << "y3(150) = " << y3(149) << ", y3(155) = " << y3(154) << endl;
        }
    }
    //----------------------------------------------------//
return 0;
}
