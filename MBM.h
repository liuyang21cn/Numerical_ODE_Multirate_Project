/*
    Yang Liu
    SMU Mathematics
    
    Multirate Base Method
 */

#ifndef MBM_DEFINED__
#define MBM_DEFINED__

// Inclusions
#include <math.h>
#include "mat.h"
#include "rhs.h"

// Forward Euler time stepper class
class MultirateBaseMethod {
    
private:
    
    // private reusable local data
    Mat *f;             // storage for ODE RHS vector
    RHSFunction *frhs;  // pointer to ODE RHS function
    
public:
    
    // constructor (sets RHS function pointer, copies y for local data)
    MultirateBaseMethod(RHSFunction &frhs_, Mat &y) {
        frhs = &frhs_;
        f = new Mat(y);
    };
    
    // destructor (frees local data)
    ~MultirateBaseMethod() { delete f; };
    
    // Evolve routine (evolves the solution via forward Euler)
    Mat Evolve(Mat tspan, double h, Mat &y, int m, int const beg, int const end);
    
    int Assemble(Mat &ynew, Mat &yslow, Mat &yfast, int const beg, int const end);
    int Split(Mat &y, Mat &yslow, Mat &yfast, int const beg, int const end);
    
    
};

#endif
