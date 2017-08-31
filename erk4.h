/* Explicit 4th-order Runge-Kutta time stepper class header file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2014  */

#ifndef ERK4_DEFINED__
#define ERK4_DEFINED__

// Inclusions
#include <math.h>
#include "mat.h"
#include "rhs.h"


// Explicit RK4 time stepper class
class ERK4Stepper {

 private:

  // private reusable local data
  RHSFunction *frhs;  // pointer to ODE RHS function

 public:

  // constructor (sets RHS function pointer)
  ERK4Stepper(RHSFunction &frhs_) { frhs = &frhs_; };

  // Evolve routine (evolves the solution)
  Mat Evolve(Mat tspan, double h, Mat &y);

  // Single step calculation
  int Step(double t, double h, Mat &z, Mat &f0, 
	   Mat &f1, Mat &f2, Mat &f3, Mat &y);

};

#endif
