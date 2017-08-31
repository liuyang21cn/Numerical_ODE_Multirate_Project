/* Explicit 4th-order Runge-Kutta time stepper class implementation file.

   D.R. Reynolds
   Math 6321 @ SMU
   Fall 2014  */

#include "mat.h"
#include "erk4.h"


// The explicit 4th-order Runge-Kutta time step evolution routine
//
// Inputs:  tspan holds the current time interval, [t0, tf]
//          h holds the desired time step size
//          y holds the initial condition, y(t0)
// Outputs: y holds the computed solution, y(tf)
//
// The return value is a row vector containing all internal 
// times at which the solution was computed,
//               [t0, t1, ..., tN]
Mat ERK4Stepper::Evolve(Mat tspan, double h, Mat &y) {

  // check for legal inputs 
  if (h <= 0.0) {
    std::cerr << "Evolve: Illegal h\n";
    return tspan;
  }
  if (tspan(1) <= tspan(0)) {
    std::cerr << "Evolve: Illegal tspan\n";
    return tspan;	  
  }
  
  // figure out how many time steps
  long int N = (tspan(1)-tspan(0))/h;
  if (tspan(1) > tspan(0)+N*h)  N++;
    
  // create ouput Mat
  Mat times(1,N+1);
  times(0) = tspan(0);

  // set temporary vectors to be reused at each step
  Mat z(y);
  Mat f0(y);
  Mat f1(y);
  Mat f2(y);
  Mat f3(y);

  // iterate over time steps
  for (int i=0; i<N; i++) {

    // last step only: update h to stop directly at final time
    if (i == N-1) 
      h = tspan(1)-times(i);

    // perform a single step of RK4 to update y
    if (Step(times(i), h, z, f0, f1, f2, f3, y) != 0) {
      std::cerr << "Evolve: Error in Step() function\n";
      return times;
    }

    // update current time, store in output array
    times(i+1) = times(i) + h;
  }

  return times;
}


// Single step of explicit 4th-order Runge-Kutta
//
// Inputs:  t holds the current time
//          h holds the current time step size
//          z, f1-f4 hold temporary vectors needed for the problem
//          y holds the current solution
// Outputs: y holds the updated solution, y(t+h)
//
// The return value is an integer indicating success/failure,
// with 0 indicating success, and nonzero failure.
int ERK4Stepper::Step(double t, double h, Mat &z, Mat &f0, 
		      Mat &f1, Mat &f2, Mat &f3, Mat &y) {

  // set Butcher tables that define the method
  double Adata[] = {0.0, 0.0, 0.0, 0.0, 
		   1.0/2.0, 0.0, 0.0, 0.0, 
		   0.0, 1.0/2.0, 0.0, 0.0, 
		   0.0, 0.0, 1.0, 0.0};
  double bdata[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
  double cdata[] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
  Mat A(4,4,Adata);  A.Trans();   // Adata is stored by column; convert to row
  Mat b(1,4,bdata);
  Mat c(1,4,cdata);

  // stage 1: set stage and compute RHS
  z = y;
  if (frhs->Evaluate(t, z, f0) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 2: set stage and compute RHS
  z = y + h*A(1,0)*f0;
  if (frhs->Evaluate(t+c(1)*h, z, f1) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 3: set stage and compute RHS
  z = y + h*A(2,0)*f0 + h*A(2,1)*f1;
  if (frhs->Evaluate(t+c(2)*h, z, f2) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // stage 4: set stage and compute RHS
  z = y + h*A(3,0)*f0 + h*A(3,1)*f1 + h*A(3,2)*f2;
  if (frhs->Evaluate(t+c(3)*h, z, f3) != 0) {
    std::cerr << "Step: Error in ODE RHS function\n";
    return 1;
  }

  // compute next step solution
  y += h*b(0)*f0 + h*b(1)*f1 + h*b(2)*f2 + h*b(3)*f3;

  // return success
  return 0;
}
