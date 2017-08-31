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
class RHS_AD: public RHSFunction {
    
public:
    
    int Evaluate(double t, Mat &y, Mat &f) {
        
        double a = 5.0;
        double d = 0.01;
        double c = 100;
        double g = 0.0;
        double pi = 3.14159265359;
        
        int Ny = y.Size();
        double dx = 2.0/(Ny-1);
        double x = 0.0;
        
        f = 0.0;
        for (int i = 2; i<=Ny-3; i++){
            
            x = -1.0+i*dx;
            g = 1000*pow(cos(pi*x/2), 200)*sin(pi*t);
            
            f(i) = ( -a*(-y(i+2)+8.0*y(i+1)-8.0*y(i-1)+y(i-2))/12.0/dx + d*(-y(i+2)+16.0*y(i+1)-30.0*y(i)+16.0*y(i-1)-y(i-2))/12.0/dx/dx - c*y(i) + g );\

        }
        
        f(0) = 0;
        f(1) = 0;
        f(Ny-1) = 0;
        f(Ny-2) = 0;
        return 0;
    }
};

// ---------------------------------------------//

// main routine
int main() {
    
    // time steps to try
    Mat h("0.001, 0.0005, 0.00025");
    
    // step for fast component
    Mat m("1, 2, 4, 16, 32, 64");
    
    //----- Advection Diffusion Problem -----//
    // initial condition
    int Ny = 401;
    double L = 1;
    
    Mat y0(Ny);
    y0 = 0.0;
    
    // matlab solution
    double y210 = 4.427312253379179;
    double y220 = 2.727795834426871;

    // set up step
    double t0 = 0.0;
    double Tf = 0.8;
    double tcur = t0;
    double dtout = 0.1;
    
    // create ODE RHS function objects
    RHS_AD rhs;
    
    // create multirate methdod objects
    MultirateBaseMethod MBM(rhs, y0);
    MPRK2Stepper        MPRK2(rhs, y0);
    MPRK3Stepper        MPRK3(rhs, y0);
    
    // temporary variables
    Mat y1(Ny), y2(Ny), y3(Ny), y4(Ny), yerr1(Ny), yerr2(Ny), yerr3(Ny), yerr4(Ny), tspan(1,2);
    double err1, maxerr1;
    Mat errs1(h.Size()), errs2(h.Size()), errs3(h.Size()), errs4(h.Size());
    
    //--------------------------------------//
    
    
    ///////// Advection-Diffusion equation /////////
    cout << "\nAdvection-Diffusion equation problem:\n";
    
    //----- Multirate Base Method -----//
    /*
     
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
                Mat tvals = MBM.Evolve(tspan, h(ih), y1, m(im), 150, 250);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
            }
            cout << "   h = " << h(ih);
            for (int ii=0; ii<401; ii++){
                cout << y1(ii) << endl;
            }
//            errs1(ih) = y1(209) - y210;
//            errs2(ih) = y1(219) - y220;
//            cout << "\t err1(210) = " << (y1(209)-y210) << ", err1(220) = " << (y1(219)-y220);
//            if (ih > 0)
//                cout << "\t  conv rate = " << (log(errs1(ih))-log(errs1(ih-1)))/(log(h(ih))-log(h(ih-1))) <<"\t, "<< (log(errs2(ih))-log(errs2(ih-1)))/(log(h(ih))-log(h(ih-1)));
//            cout << endl;
        }
    }
    //--------------------------------------------//
    */
    //----- Constantinescu and Sandu RK2a Method -----//
   /* cout << "\n 2. Constantinescu and Sandu RK2a Method:\n";
    
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
                Mat tvals = MPRK2.Evolve(tspan, h(ih), y2, m(im), 150, 250);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
            }
            cout << "   h = " << h(ih);
            for (int ii=0; ii<401; ii++){
                cout << y2(ii) << endl;
            }
//            cout << "   h = " << h(ih);
//            errs3(ih) = y2(209) - y210;
//            errs4(ih) = y2(219) - y220;
//            cout << "\t err3(210) = " << errs3(ih) << ", err4(220) = " << errs4(ih);
//            if (ih > 0)
//                cout << "\t  conv rate = " << (log(errs3(ih))-log(errs3(ih-1)))/(log(h(ih))-log(h(ih-1))) <<"\t, "<< (log(errs4(ih))-log(errs4(ih-1)))/(log(h(ih))-log(h(ih-1)));
//            cout << endl;
        }
    }
    //----------------------------------------------------//
    
    //----- Constantinescu and Sandu RK43 Method -----//
   */ cout << "\n 3. Constantinescu and Sandu RK43 Method:\n";
    
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
                Mat tvals = MPRK3.Evolve(tspan, h(ih), y3, m(im), 150, 250);
                tcur = tvals(tvals.Size()-1);   // last entry in tvals
                
            }
            cout << "   h = " << h(ih);
            errs1(ih) = y3(209) - y210;
            errs2(ih) = y3(219) - y220;
            cout << "\t err1(210) = " << (y3(209)-y210) << ", err2(220) = " << (y3(219)-y220);
            if (ih > 0)
                cout << "\t  conv rate = " << (log(errs1(ih))-log(errs1(ih-1)))/(log(h(ih))-log(h(ih-1))) << "\t, " << (log(errs2(ih))-log(errs2(ih-1)))/(log(h(ih))-log(h(ih-1)));
            cout << endl;
        }
    }
    //----------------------------------------------------//*/
    return 0;
}
