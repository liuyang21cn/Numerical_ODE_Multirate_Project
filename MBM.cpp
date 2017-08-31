/*
 
 Yang Liu
 SMU Mathematics
 
 Multirate Base Method
 */

#include "mat.h"
#include "MBM.h"
#include "rhs.h"
using namespace std;

// Assemble function is to put slow and fast parts into y
int MultirateBaseMethod::Assemble(Mat &ynew, Mat &yslow, Mat &yfast, int const beg, int const end){
    
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
int MultirateBaseMethod::Split(Mat &y, Mat &yslow, Mat &yfast, int const beg, int const end){
    
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


Mat MultirateBaseMethod::Evolve(Mat tspan, double h1, Mat &y, int m, int beg, int end) {
    
    // check for legal inputs
    if (h1 <= 0.0) {
        std::cerr << "Multirate Base Method: Illegal h\n";
        return tspan;
    }
    if (tspan(1) <= tspan(0)) {
        std::cerr << "Multirate Base Method: Illegal tspan\n";
        return tspan;
    }
    
    // figure out how many time steps
    long int N = (tspan(1)-tspan(0))/h1;
    if (tspan(1) > tspan(0)+N*h1)  N++;
    
    int N2 = m;
    long double h2 = h1/m;
    
    // create ouput Mat for fast and slow parts
    Mat times1(1, N+1);
    times1(0) = tspan(0);
    
    Mat times2(1, N2+1);
    times2(0) = tspan(0);
    
    // beg and end use C++ format indices, as 0 is the first index of array
    long int nfast = end-beg+1;
    long int nslow = y.Size()-nfast;
    
    Mat yfast(nfast);
    Mat yslow(nslow);

    Mat yslow1(nslow), yslow2(nslow);
    Mat fslow(nslow), ffast(nfast);
    Mat temps(nslow), tempf(nfast), ytemp(y.Size());

    
    // iterate over time steps
   for (int i=0; i<N; i++) {
    
       // update yslow and yfast
//    cout << "begin split" << endl;
        Split(y, yslow, yfast, beg, end);

       // update y(n)
        yslow1 = yslow;

       // last step only: update h to stop directly at final time

       if (i == N-1)   h1 = tspan(1)-times1(i);

       // update time spacing for fast parts
       h2 = h1/N2;

    // compute ODE RHS for slow parts
        if (frhs->Evaluate(times1(i), y, *f) != 0) {
            std::cerr << "Multirate Base Method: Error in 1st ODE RHS function\n";
            return times1;
        }

       // update solution with forward Euler step
       Split(*f, fslow, ffast, beg, end);
        yslow += h1*(fslow);

       // update y(n+1)
      yslow2 = yslow;

       // update current time, store in output array
       times1(i+1) = times1(i) + h1;

       // cout << "times1(i) = " << times1(i) << endl;

       // update start time of fast parts be the current time for slow parts
        times2(0) = times1(i);

    // evaluate fast component
      for (int i2=0; i2<N2; i2++) {

          // last step only: update h to stop directly at final time
           if (i2 == N2-1)
               h2 = times1(i+1)-times2(i2);

          // update y temp
            temps = ( m - (i2+1.0) + 1.0)/m*yslow1 + ((i2+1.0) - 1.0)/m*yslow2;
            tempf = yfast;
          
            Assemble(ytemp, temps, tempf, beg, end);
          
            // compute ODE RHS for fast parts
            if (frhs->Evaluate(times2(i2), ytemp, *f) != 0) {
                std::cerr << "Multirate Base Method: Error in 2nd ODE RHS function\n";
                return times2;
            }

          // update solution with forward Euler step
            Split(*f, fslow, ffast, beg, end);
            yfast += h2*(ffast);

          // update current time, store in output array
           times2(i2+1) = times2(i2) + h2;
           //cout << "times2(i) = "  << times2(i2) <<endl;
        }

       // update result y
        Assemble(y, yslow, yfast, beg, end);

   }
    
    return times1;
}

