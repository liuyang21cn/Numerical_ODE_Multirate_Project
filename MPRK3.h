/*
  
  Yang Liu
  SMU Mathematics

  
*/


#ifndef MPRK3_DEFINED__
#define MPRK3_DEFINED__

// Inclusions
#include <math.h>
#include "rhs.h"

class MPRK3Stepper {
    
private:
    
    // private reusable local data
    Mat *f;             // storage for ODE RHS vector
    RHSFunction *frhs;  // pointer to ODE RHS function
    
public:
    
    // constructor (sets RHS function pointer)
    MPRK3Stepper(RHSFunction &frhs_, Mat &y) {
        frhs = &frhs_;
        f = new Mat(y);
    };
    
    // destructor (frees local data)
    ~MPRK3Stepper() { delete f;};
    
    // Evolve routine (evolves the solution)
    Mat Evolve(Mat tspan, double h, Mat &y, int m, int const beg, int const end);
    
    // Assembel routine to pass information y
    int Assemble(Mat &ynew, Mat &yslow, Mat &yfast, int const beg, int const end);
    
    // Split routine to pass information to yfast & yslow
    int Split(Mat &y, Mat &yslow, Mat &yfast, int const beg, int const end);
    
};

#endif

