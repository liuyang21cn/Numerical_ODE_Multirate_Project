/* ODE RHS and Jacobian function abstract base class definitions.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2014  */

#ifndef ODE_RHS_DEFINED__
#define ODE_RHS_DEFINED__

// Inclusions
#include "mat.h"


// Declare abstract base classes for ODE RHS and its Jacobian, to 
// define what the backward Euler solver expects from each.

//   ODE RHS function abstract base class; derived classes 
//   must at least implement the Evaluate() routine
class RHSFunction {
 public: 
  virtual int Evaluate(double t, Mat &y, Mat &f) = 0;
};

//   ODE RHS Jacobian function abstract base class; derived 
//   classes must at least implement the Evaluate() routine
class RHSJacobian {
 public: 
  virtual int Evaluate(double t, Mat &y, Mat &J) = 0;
};

#endif
