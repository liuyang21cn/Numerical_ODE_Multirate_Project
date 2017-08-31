/* Daniel R. Reynolds
   SMU Mathematics
   13 August 2014 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <string.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
using namespace std;

// singularity tolerance
#define STOL 1.e-15


// This file implements the operations defined in the Mat class


///// General class routines /////

// matrix constructor (initializes values to 0.0)
Mat::Mat(long int m, long int n) {
  nrows = m;
  ncols = n;
  data = new double[n*m];
  own_data = true;
  for (long int i=0; i<n*m; i++)  data[i] = 0.0;
}

// column vector constructor (initializes values to 0.0)
Mat::Mat(long int m) {
  nrows = m;
  ncols = 1;
  data = new double[m];
  own_data = true;
  for (long int i=0; i<m; i++)  data[i] = 0.0;
}

// string utility routines we'll need
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

// string-based matrix constructor
Mat::Mat(string mat_spec) {
  // initialize empty matrix
  nrows = 0;
  ncols = 0;
  data = NULL;
  own_data = false;

  // parse string to determine matrix rows
  vector<string> rows = split(mat_spec, ';');
  long int m = rows.size();
  if (m < 1) {
    cerr << "string Mat constructor error: empty string!\n";
    return;
  }

  // verify that all rows have the same number of columns
  vector<string> row = split(rows[0], ',');
  long int n = row.size();
  for (long int i=1; i<nrows; i++) {
    row = split(rows[i], ',');
    if (row.size() != n) {
      cerr << "string Mat constructor error: all rows must have the same number of columns!\n";
      return;
    }
  }
   
  // allocate matrix data
  data = new double[n*m];

  // fill Mat structure
  own_data = true;
  nrows = m;
  ncols = n;
  for (long int i=0; i<nrows; i++) {
    row = split(rows[i], ',');
    for (long int j=0; j<ncols; j++) {
      stringstream ss(row[j]);
      ss >> data[IDX(i,j,nrows)];
    }
  }
}

// empty constructor (sets data to point to existing array)
Mat::Mat(long int m, long int n, double *data_array) {
  nrows = m;
  ncols = n;
  data = data_array;
  own_data = false;
}

// copy constructor
Mat::Mat(const Mat& A) {
  nrows = A.Rows();
  ncols = A.Cols();
  data = new double[ncols*nrows];
  own_data = true;
  for (long int i=0; i<nrows*ncols; i++)  data[i] = A.data[i];
}

// write myself to stdout
int Mat::Write() const {

  // return with failure if data array isn't allocated
  if (data == NULL) {
    cerr << "Mat:Write error, empty data array\n";
    return 1;
  }

  // print data to screen 
  for (long int i=0; i<nrows; i++) {
    for (long int j=0; j<ncols; j++)
      printf("  %.17g",data[IDX(i,j,nrows)]);
    printf("\n");
  }
  return 0;
}


// write myself to a file
int Mat::Write(const char *outfile) const {

  // return with failure if data array isn't allocated
  if (data == NULL) {
    cerr << "Mat:Write error, empty data array\n";
    return 1;
  }

  // return with failure if 'outfile' is empty
  if (strlen(outfile) < 1) {
    cerr << "Mat::Write error, empty outfile\n";
    return 1;
  }

  // open output file
  FILE *fptr = NULL;
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    cerr << "Mat::Write error, unable to open " << outfile << " for writing\n";
    return 1;
  }

  // print data to file
  for (long int i=0; i<nrows; i++) {
    for (long int j=0; j<ncols; j++)
      fprintf(fptr, "  %.16g",data[IDX(i,j,nrows)]);
    fprintf(fptr, "\n");
  }

  // close output file and return
  fclose(fptr);
  return 0;
}


// streaming output operator
ostream& operator<<(ostream &out, Mat &data) {

  // print data to "out" stream and return
  for (long int i=0; i<data.Rows(); i++) {
    for (long int j=0; j<data.Cols(); j++)
      out << "  " << data(i,j);
    out << endl;
  }
  return out;
}


// access my data -- access value by location in matrix
// gets reference to value so that it can be changed
double& Mat::operator()(long int i, long int j) {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat() error: empty data array\n";
    return data[0];
  }

  // check that (i,j) is a valid entry
  if (i<0 || i>=nrows || j<0 || j>=ncols) {
    cerr << "Mat() error: (" << i << "," << j << ") outside matrix bounds\n";
    cerr << "  returning first entry to avoid segmentation fault\n";
    return data[0];
  }
  return data[IDX(i,j,nrows)];
}


// access my data -- access value by location in data array,
// gets reference to value so that it can be changed
double& Mat::operator()(long int i){

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat() error: empty data array\n";
    return data[0];
  }

  // check that (i) is a valid entry
  if (i<0 || i>=nrows*ncols) {
    cerr << "Mat() error: (" << i << ") outside vector bounds\n";
    cerr << "  returning first entry to avoid segmentation fault\n";
    return data[0];
  }
  return data[i];
}



///// Arithmetic operations defined on a given Mat /////


// C = A*x + C
int Mat::Axpy(const Mat &A, const Mat &x) {

  // check that array sizes are acceptable
  if ((nrows != A.nrows) || (A.ncols != x.nrows) ||
      (ncols != 1) || (x.ncols != 1)) {
    cerr << "Mat::Axpy error, illegal matrix dimensions\n";
    cerr << "  A is " << A.nrows << " x " << A.ncols << ",  x is " << x.nrows 
	 << " x " << x.ncols << ", result is " << nrows << " x " 
	 << ncols << endl;
    return 1;
  }
  
  // do operation
  for (long int j=0; j<A.ncols; j++)
    for (long int i=0; i<nrows; i++) 
      data[i] += A.data[IDX(i,j,nrows)] * x.data[j];
  
  // return with success
  return 0;
}


// C = a*C + b*B
int Mat::LinearSum(double a, double b, const Mat &B) {

  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    cerr << "Mat::LinearSum error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ", Mat 2 is " << B.nrows << " x " << B.ncols << endl;
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || B.data == NULL) {
    cerr << "Mat::LinearSum error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = a*data[i] + b*B.data[i];
  
  // return success
  return 0;
}


// C = a*A + b*B
int Mat::LinearSum(double a, const Mat &A, double b, const Mat &B) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols || 
      B.nrows != nrows || B.ncols != ncols) {
    cerr << "Mat::LinearSum error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << A.nrows << " x " << A.ncols 
	 << ",  Mat 3 is " << B.nrows << " x " << B.ncols << endl;
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL || B.data == NULL) {
    cerr << "Mat::LinearSum error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = a*A.data[i] + b*B.data[i];
  
  // return success
  return 0;
}


//   C = C+a  (adds scalar a to my data)
int Mat::Add(double a) {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat::Add error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] += a;
  
  // return success
  return 0;
}


//   C = C.*A (component-wise multiply of my data by A)
int Mat::Mul(const Mat &A) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Mat::Mul error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << A.nrows << " x " << A.ncols << endl;
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL) {
    cerr << "Mat::Mul error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] *= A.data[i];
  
  // return success
  return 0;
}


//   C = C./A (component-wise division of my data by A)
int Mat::Div(const Mat &A) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Mat::Div error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << A.nrows << " x " << A.ncols << endl;
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL) {
    cerr << "Mat::Div error: empty data array\n";
    return 1;
  }

  // check that A is nowhere zero, otherwise return with error
  for (long int i=0; i<nrows*ncols; i++)
    if (A.data[i] == 0.0) 
      return 1;

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] /= A.data[i];
  
  // return success
  return 0;
}


//   C = a*C  (scales my data by scalar a)
int Mat::Mul(double a) {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat::Mul error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  data[i] *= a;
  
  // return success
  return 0;
}


//   C = A  (copies A into C)
int Mat::Copy(const Mat &A) {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Mat::Copy error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << A.nrows << " x " << A.ncols << endl;
    return 1;
  }
  
  // check that data is not NULL
  if (data == NULL || A.data == NULL) {
    cerr << "Mat::Copy error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = A.data[i];
  
  // return success
  return 0;
}


//   C(is:ie,js:je) = A  (copies A into a submatrix of C)
int Mat::Copy(const Mat &A, long int is, long int ie, 
	      long int js, long int je) {

  // check that array sizes match
  if (A.nrows != (ie-is+1) || A.ncols != (je-js+1)) {
    cerr << "Mat::Copy error, matrix size mismatch\n";
    cerr << "  supplied Mat is " << A.nrows << " x " << A.ncols 
	 << ", but requested submatrix is " << ie-is+1 << " x " 
	 << je-js+1 << endl;
    return 1;
  }
  // check for valid submatrix
  if (is < 0 || is >= nrows) {
    cerr << "Mat::Copy error, requested submatrix does not exist\n";
    cerr << "  illegal is = " << is << " (matrix has " << nrows << " rows)\n";
    return 1;
  }
  if (ie < 0 || ie >= nrows) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal ie = " << ie << " (matrix has " << nrows << " rows)\n";
    return 1;
  }
  if (js < 0 || js >= ncols) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal js = " << js << " (matrix has " << ncols << " columns)\n";
    return 1;
  }
  if (je < 0 || je >= ncols) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal je = " << je << " (matrix has " << ncols << " columns)\n";
    return 1;
  }
  if (ie < is || je < js) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  upper index is below lower index: is = " << is << ", ie = " 
	 << ie << ", js = " << js << ", je = " << je << endl;
    return 1;
  }
  
  // check that data arrays are not NULL
  if (data == NULL || A.data == NULL) {
    cerr << "Mat::Copy error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int j=0; j<A.ncols; j++)  
    for (long int i=0; i<A.nrows; i++)  
      data[IDX(i+is,j+js,nrows)] = A.data[IDX(i,j,A.nrows)];
  
  // return success
  return 0;
}


//   C = a  (sets all entries of C to the scalar a)
int Mat::Const(double a) {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat::Const error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = a;
  
  // return success
  return 0;
}


// C = C.^p (component-wise exponentiation of my data to the power p)
int Mat::Power(double p) {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat::Power error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = pow(data[i], p);
  
  // return success
  return 0;
}


// Cij = |Cij| (component-wise absolute value of my data)
int Mat::Abs() {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat::Abs error: empty data array\n";
    return 1;
  }

  // perform operation
  for (long int i=0; i<nrows*ncols; i++)  
    data[i] = fabs(data[i]);
  
  // return success
  return 0;
}


// Cij = Cji
int Mat::Trans() {

  // check that data is not NULL
  if (data == NULL) {
    cerr << "Mat::Trans error: empty data array\n";
    return 1;
  }

  // perform operation in place if matrix is square
  if (nrows == ncols) {
    double tmp;
    for (long int i=0; i<nrows; i++)
      for (long int j=0; j<i; j++) {
	tmp = data[IDX(i,j,nrows)];
	data[IDX(i,j,nrows)] = data[IDX(j,i,nrows)];
	data[IDX(j,i,nrows)] = tmp;
      }

  // otherwise we need a new data array to ensure a clean transpose
  } else {

    // create temporary data array, and copy transposed data over 
    double *newdata = new double[nrows*ncols];
    for (long int i=0; i<nrows; i++)
      for (long int j=0; j<ncols; j++)
	newdata[IDX(j,i,ncols)] = data[IDX(i,j,nrows)];

    // copy newdata values into existing data array
    for (long int i=0; i<nrows*ncols; i++)
      data[i] = newdata[i];

    // delete temporary data array
    delete[] newdata;
    
    // swap matrix dimensions
    int tmp = nrows;
    nrows = ncols;
    ncols = tmp;
  }
  
  // return success
  return 0;
}


// fill in an existing vector with linearly spaced data
int Mat::Linspace(double a, double b) {

  // ensure that this is a row or column vector
  if (nrows != 1 && ncols != 1) {
    cerr << "Mat::Linspace error, matrix must be a row/column vector\n";
    cerr << "  dimensions are " << nrows << " x " << ncols << endl;
    return 1;
  }

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = a + (b-a)/(nrows*ncols-1)*i;
  return 0;
}


// fill in an existing vector with logarithmically spaced data
int Mat::Logspace(double a, double b) {

  // ensure that this is a row or column vector
  if (nrows != 1 && ncols != 1) {
    cerr << "Mat::Logspace error, matrix must be a row/column vector\n";
    cerr << "  dimensions are " << nrows << " x " << ncols << endl;
    return 1;
  }

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = pow(10.0, a + (b-a)/(nrows*ncols-1)*i);
  return 0;
}


// fill in a vector with uniformly-distributed random numbers in [0,1]
int Mat::Random() {

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = random() / (pow(2.0,31.0) - 1.0);
  return 0;
}


// fill in a matrix as the identity
int Mat::Eye() {

  // ensure that this is a square matrix
  if (nrows != ncols) {
    cerr << "Mat::Eye error, matrix must be square\n";
    cerr << "  dimensions are " << nrows << " x " << ncols << endl;
    return 1;
  }

  // fill in entries and return
  for (long int i=0; i<nrows*ncols; i++)
    data[i] = 0.0;
  for (long int i=0; i<nrows; i++)
    data[IDX(i,i,nrows)] = 1.0;
  return 0;
}



// computes the inverse of a nonsingular matrix 
int Mat::Inverse() {

  // check that matrix sizes match
  if (nrows != ncols) {
    cerr << "Inverse error, non-square matrix\n";
    cerr << "  Mat is " << nrows << " x " << ncols << endl;
    return 1;
  }

  // create two temporary matrices for operation
  Mat B(nrows,ncols);
  Mat A(nrows,ncols);

  // fill RHS matrix with I, and A matrix with myself
  B.Eye();
  A.Copy(*this);

  // call existing Solve routine for computations
  if (Solve(A, *this, B) != 0)
    return 1;

  // return success
  return 0;
}


///// New matrix creation operations (A is calling the operation, C is the result) /////



// C = A+B
Mat& Mat::operator+(const Mat &B) {

  // create a pointer to the result
  Mat *C = NULL;

  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    cerr << "Mat::operator+ error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << B.nrows << " x " << B.ncols << endl;
    C = new Mat(0,0);
    return *C;
  }

  // create new Mat for output, and do operation
  C = new Mat(nrows,ncols);
  for (long int i=0; i<nrows*ncols; i++)  
    C->data[i] = data[i] + B.data[i];

  // return result
  return *C;
}


// C = A-B
Mat& Mat::operator-(const Mat &B) {

  // create a pointer to the result
  Mat *C = NULL;

  // check that array sizes match
  if (B.nrows != nrows || B.ncols != ncols) {
    cerr << "Mat::operator- error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << B.nrows << " x " << B.ncols << endl;
    C = new Mat(0,0);
    return *C;
  }

  // create new Mat for output, and do operation
  C = new Mat(nrows,ncols);
  for (long int i=0; i<nrows*ncols; i++)  
    C->data[i] = data[i] - B.data[i];

  // return result
  return *C;
}


// C = A*B
Mat& Mat::operator*(const Mat &B) {

  // create a pointer to the result
  Mat *C = NULL;

  // determine if either matrix is a scalar
  bool A_scalar = ((nrows==1) && (ncols==1));
  bool B_scalar = ((B.nrows==1) && (B.ncols==1));

  // scalar-times-matrix
  if (A_scalar) {

    // create new Mat for output, and do operation
    C = new Mat(B.nrows,B.ncols);
    for (long int i=0; i<nrows*ncols; i++)  
      C->data[i] = data[0] * B.data[i];
    
  // scalar-times-matrix
  } else if (B_scalar) {
    
    // create new Mat for output, and do operation
    C = new Mat(nrows,ncols);
    for (long int i=0; i<nrows*ncols; i++)  
      C->data[i] = data[i] * B.data[0];
        
  // normal matrix product
  } else {

    // check that array sizes are acceptable
    if (B.nrows != ncols) {
      cerr << "Mat::operator* error, inner dimension mismatch\n";
      cerr << "  Mat 1 is " << nrows << " x " << ncols 
	   << ",  Mat 2 is " << B.nrows << " x " << B.ncols << endl;
      C = new Mat(0,0);
      return *C;
    }
    
    // create new Mat for output, and do operation
    C = new Mat(nrows,B.ncols);
    for (long int i=0; i<nrows; i++) 
      for (long int j=0; j<B.ncols; j++)
	for (long int k=0; k<ncols; k++)
	  C->data[IDX(i,j,nrows)] += data[IDX(i,k,nrows)] * B.data[IDX(k,j,B.nrows)];

  }

  // return result
  return *C;
}


// C = A*b
Mat& Mat::operator*(const double b) {

  // create a pointer to the result
  Mat *C = new Mat(nrows,ncols);

  for (long int i=0; i<nrows*ncols; i++)  
    C->data[i] = data[i] * b;
    
  // return result
  return *C;
}


// C = a*B
Mat& operator*(const double a, const Mat &B) {

  // create a pointer to the result
  Mat *C = new Mat(B.Rows(),B.Cols());
  C->LinearSum(0.0, a, B);
    
  // return result
  return *C;
}


// C = A^T
Mat Mat::T() {

  // create new Mat for output, and do operation
  Mat *C = new Mat(ncols,nrows);
  for (long int j=0; j<ncols; j++)  
    for (long int i=0; i<nrows; i++)  
      C->data[IDX(j,i,ncols)] = data[IDX(i,j,nrows)];

  // return result
  return *C;
}


// column copy routine (creates a column vector from one column of an existing Mat)
Mat Mat::Column(long int j) {

  // create a pointer to the result
  Mat *C = NULL;

  // check that requested column exists
  if (j < 0 || j > ncols) {
    cerr << "Mat::GetColumn error, requested column does not exist\n";
    cerr << "  j is " << j << ", but matrix only has " << ncols << " columns\n";
    C = new Mat(0,0);
    return *C;
  }

  // create new matrix of desired size
  C = new Mat(nrows, 1);

  // copy requested data
  for (long int i=0; i<nrows; i++) 
    C->data[i] = data[IDX(i,j,nrows)];

  // return object
  return *C;
}


// row copy routine (creates a row vector from one row of an existing Mat)
Mat Mat::Row(long int i) {

  // create a pointer to the result
  Mat *C = NULL;

  // check that requested row exists
  if (i < 0 || i > nrows) {
    cerr << "Mat::GetRow error, requested row does not exist\n";
    cerr << "  i is " << i << ", but matrix only has " << nrows << " rows\n";
    C = new Mat(0,0);
    return *C;
  }

  // create new matrix of desired size
  C = new Mat(1, ncols);

  // copy requested data
  for (long int j=0; j<ncols; j++) 
    C->data[j] = data[IDX(i,j,nrows)];

  // return object
  return *C;
}


// submatrix copy routine (creates a matrix from a portion of an existing Mat)
Mat Mat::operator()(long int is, long int ie, long int js, long int je) {

  // create a pointer to the result
  Mat *C = NULL;

  // check that requested submatrix exists
  if (is < 0 || is >= nrows) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal is = " << is << " (matrix has " << nrows << " rows)\n";
    C = new Mat(0,0);
    return *C;
  }
  if (ie < 0 || ie >= nrows) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal ie = " << ie << " (matrix has " << nrows << " rows)\n";
    C = new Mat(0,0);
    return *C;
  }
  if (js < 0 || js >= ncols) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal js = " << js << " (matrix has " << ncols << " columns)\n";
    C = new Mat(0,0);
    return *C;
  }
  if (je < 0 || je >= ncols) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  illegal je = " << je << " (matrix has " << ncols << " columns)\n";
    C = new Mat(0,0);
    return *C;
  }
  if (ie < is || je < js) {
    cerr << "Mat::submatrix copy error, requested submatrix does not exist\n";
    cerr << "  upper index is below lower index: is = " << is << ", ie = " 
	 << ie << ", js = " << js << ", je = " << je << endl;
    C = new Mat(0,0);
    return *C;
  }

  // create new matrix of desired size
  C = new Mat(ie-is+1, je-js+1);

  // copy requested data
  for (long int j=js; j<=je; j++) 
    for (long int i=is; i<=ie; i++) 
      C->data[IDX(i-is,j-js,ie-is+1)] = data[IDX(i,j,nrows)];

  // return object
  return *C;
}


// column accessor routine (does not alocate its own data, only returns a 
// Mat object whose data array points to an existing Mat column).
Mat* Mat::AccessColumn(long int j) {

  // create a pointer to the result
  Mat *C = NULL;

  // check that requested column exists
  if (j < 0 || j > ncols) {
    cerr << "Mat::AccessColumn error, requested column does not exist\n";
    cerr << "  j is " << j << ", but matrix only has " << ncols << " columns\n";
    C = new Mat(0,0);
  }

  // create new matrix, pointing to column j of this matrix
  C = new Mat(nrows, 1, &(data[IDX(0,j,nrows)]));

  // return result (the pointer to the object)
  return C;
}



///// Scalar output operators on matrices /////



// vector 2-norm, square root of the matrix Frobenius norm
double Mat::Norm() {

  // check that my data is allocated
  if (data == NULL) {
    cerr << "Mat::Norm error, data not allocated\n";
    return -1.0;
  }
  
  // perform operation
  double sum=0.0;
  for (long int i=0; i<nrows*ncols; i++)  sum += data[i]*data[i];
  return sqrt(sum);
}


// vector infinity norm, largest absolute value entry of a matrix
double Mat::MaxNorm() {

  // check that my data is allocated
  if (data == NULL) {
    cerr << "Mat::MaxNorm error, data not allocated\n";
    return -1.0;
  }
  
  // perform operation
  double mx=0.0;
  for (long int i=0; i<nrows*ncols; i++)  
    mx = (mx > fabs(data[i])) ? mx : fabs(data[i]);
  return mx;
}


// row vector one norm, column vector infinity norm, matrix infinity norm
double Mat::InfNorm() {

  // check that my data is allocated
  if (data == NULL) {
    cerr << "Mat::InfNorm error, data not allocated\n";
    return -1.0;
  }
  
  // perform operation
  double mx=0.0;
  for (long int i=0; i<nrows; i++) {
    double sum=0.0;
    for (long int j=0; j<ncols; j++) 
      sum += fabs(data[IDX(i,j,nrows)]);
    mx = (mx > sum) ? mx : sum;
  }
  return mx;
}
 

// row vector infinity norm, column vector one norm, matrix one norm
double Mat::OneNorm() {

  // check that my data is allocated
  if (data == NULL) {
    cerr << "Mat::OneNorm error, data not allocated\n";
    return -1.0;
  }
  
  // perform operation
  double mx=0.0;
  for (long int j=0; j<ncols; j++) {
    double sum=0.0;
    for (long int i=0; i<nrows; i++)
      sum += fabs(data[IDX(i,j,nrows)]);
    mx = (mx > sum) ? mx : sum;
  }
  return mx;
}


// min Cij
double Mat::Min() {

  // check that my data is allocated
  if (data == NULL) {
    cerr << "Mat::Min error, data not allocated\n";
    return -1.0;
  }
  
  // perform operation
  double mn=data[0];
  for (long int i=0; i<nrows*ncols; i++)  
    mn = (mn < data[i]) ? mn : data[i];
  return mn;
}


// max Cij
double Mat::Max() {

  // check that my data is allocated
  if (data == NULL) {
    cerr << "Mat::Max error, data not allocated\n";
    return -1.0;
  }
  
  // perform operation
  double mx=data[0];
  for (long int i=0; i<nrows*ncols; i++)  
    mx = (mx > data[i]) ? mx : data[i];
  return mx;
}


// equivalence-checking operator
bool Mat::operator==(Mat &A)  {

  // check that array sizes match
  if (A.nrows != nrows || A.ncols != ncols) {
    cerr << "Mat::operator== error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << nrows << " x " << ncols 
	 << ",  Mat 2 is " << A.nrows << " x " << A.ncols << endl;
    return false;
  }

  // perform operation and return
  bool equal=true;
  for (long int j=0; j<ncols; j++)  
    for (long int i=0; i<nrows; i++)
      equal &= (A(i,j) == data[IDX(i,j,nrows)]);
  return equal;
}



///////// independent Mat routines /////////



// create a new row vector of linearly spaced data
Mat Linspace(double a, double b, long int n) {

  // create vector, call in-place Linspace routine, and return
  Mat *x = new Mat(1,n);
  x->Linspace(a,b);
  return *x;
}


// create a new row vector of logarithmically spaced data
Mat Logspace(double a, double b, long int n) {

  // create vector, call in-place Logspace routine, and return
  Mat *x = new Mat(1,n);
  x->Logspace(a,b);
  return *x;
}


// create a new matrix of uniformly-distributed random data in [0,1]
Mat Random(long int m, long int n) {

  // create matrix, call in-place Random routine, and return
  Mat *R = new Mat(m,n);
  R->Random();
  return *R;
}


// create a new n by n identity matrix
Mat Eye(long int n) {

  // create matrix, call in-place Eye routine, and return
  Mat *I = new Mat(n,n);
  I->Eye();
  return *I;
}


// creates a matrix from a specified input file
Mat Read(const char *infile) {
  
  // determine matrix size
  long int nrows=0, ncols=0;
  ifstream ifs;
  string line;
  ifs.open(infile);
  while (getline(ifs, line)) {
    istringstream iss(line);   // convert line to stringstream
    float value;               // determine the number of columns on this row
    long int n=0;
    while (iss >> value)  n++;
    if ((n > 0) && (nrows == 0))  // first row, set ncols
      ncols = n;
    if ((n > 0) && (n != ncols)) {  // later row, with bad number of columns
      cerr << "Read() error, not all rows in file " << infile 
	   << " have the same number of cols, "
	   << n << " != " << ncols << endl;
      Mat *Afail = new Mat(0,0);
      return *Afail;
    }
    if (n > 0) nrows++;          // legal row, increment counter
  }
  ifs.close();

  // create matrix of desired size
  Mat *A = new Mat(nrows,ncols);

  // load matrix based on data from file
  ifs.open(infile);   // reopen input file
  for (long int i=0; i<nrows; i++) {
    getline( ifs, line ); 
    istringstream iss(line);   // convert line to stringstream
    for (long int j=0; j<ncols; j++) 
      iss >> (*A)(i,j);
  }
  ifs.close();

  // return result
  return *A;
}


// compute the dot-product of two compatible vectors x and y
double Dot(Mat &x, Mat &y) {

  // check that array sizes match
  if (y.Rows() != x.Rows() || y.Cols() != x.Cols()) {
    cerr << "Dot error, matrix size mismatch\n";
    cerr << "  Mat 1 is " << x.Rows() << " x " << x.Cols() 
	 << ",  Mat 2 is " << y.Rows() << " x " << y.Cols() << endl;
    return 0.0;
  }
  
  // perform operation and return
  double sum=0.0;
  for (long int j=0; j<x.Cols(); j++)  
    for (long int i=0; i<x.Rows(); i++)  
      sum += x(i,j)*y(i,j);
  return sum;
}


// performs backwards substitution on the linear system U*x = b, filling in the input Mat x
int BackSub(Mat &Umat, Mat &Xmat, Mat &Bmat) {

  // check that matrix sizes match
  if (Umat.Rows() != Bmat.Rows() || Umat.Rows() != Umat.Cols() || 
      Bmat.Cols() != Xmat.Cols() || Xmat.Rows() != Umat.Rows()) {
    cerr << "BackSub error, illegal matrix/vector dimensions\n";
    cerr << "  Mat is " << Umat.Rows() << " x " << Umat.Cols() 
	 << ",  rhs is " << Bmat.Rows() << " x " << Bmat.Cols()
	 << ",  solution is " << Xmat.Rows() << " x " << Xmat.Cols() << endl;
    return 1;
  }
  
  // get the matrix size 
  long int n = Umat.Rows();
  
  // access the data arrays
  double *U = Umat.get_data();
  double *X = Xmat.get_data();

  // copy b into x
  Xmat = Bmat;

  // analyze matrix for typical nonzero magnitude
  double Umax = Umat.MaxNorm();

  // perform column-oriented Backwards Substitution algorithm
  for (long int j=n-1; j>=0; j--) {

    // check for nonzero matrix diagonal
    if (fabs(U[IDX(j,j,n)]) < STOL*Umax) {
      cerr << "BackSub error: numerically singular matrix!\n";
      return 1;
    }

    // solve for this row of solution
    for (long int k=0; k<Xmat.Cols(); k++) 
      X[IDX(j,k,n)] /= U[IDX(j,j,n)];

    // update all remaining rhs
    for (long int i=0; i<j; i++)
      for (long int k=0; k<Xmat.Cols(); k++)
	X[IDX(i,k,n)] -= U[IDX(i,j,n)]*X[IDX(j,k,n)];

  }

  // return success
  return 0;
}


// performs backwards substitution on the linear system U*x = b, returning x as a new Mat; 
// leaves U and b untouched
Mat BackSub(Mat &U, Mat &B) {

  // check that matrix sizes match
  if (U.Rows() != B.Rows() || U.Rows() != U.Cols()) {
    cerr << "BackSub error, illegal matrix/vector dimensions\n";
    cerr << "  Mat is " << U.Rows() << " x " << U.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols() << endl;
    Mat *X = new Mat(0,0);
    return *X;
  }
  
  // create new Mat for output
  Mat *X = new Mat(U.Rows(),B.Cols());

  // call existing BackSub routine for computations
  if (BackSub(U, *X, B) != 0)
    cerr << "BackSub Warning: error in BackSub call\n";

  // return result
  return *X;
}


// performs forwards substitution on the linear system L*x = b, filling in the input Mat x
int FwdSub(Mat &Lmat, Mat &Xmat, Mat &Bmat) {

  // check that matrix sizes match
  if (Lmat.Rows() != Bmat.Rows() || Lmat.Rows() != Lmat.Cols() || 
      Bmat.Cols() != Xmat.Cols() || Xmat.Rows() != Lmat.Rows()) {
    cerr << "FwdSub error, illegal matrix/vector dimensions\n";
    cerr << "  Mat is " << Lmat.Rows() << " x " << Lmat.Cols() 
	 << ",  rhs is " << Bmat.Rows() << " x " << Bmat.Cols()
	 << ",  solution is " << Bmat.Rows() << " x " << Bmat.Cols() << endl;
    return 1;
  }
  
  // get the matrix size 
  long int n = Lmat.Rows();
  
  // access the data arrays
  double *L = Lmat.get_data();
  double *X = Xmat.get_data();

  // copy B into X
  Xmat = Bmat;

  // analyze matrix for typical nonzero magnitude
  double Lmax = Lmat.MaxNorm();

  // perform column-oriented Forwards Substitution algorithm
  for (long int j=0; j<n; j++) {

    // check for nonzero matrix diagonal
    if (fabs(L[IDX(j,j,n)]) < STOL*Lmax) {
      cerr << "FwdSub error: singular matrix!\n";
      return 1;
    }

    // solve for this row of solution
    for (long int k=0; k<Xmat.Cols(); k++)
      X[IDX(j,k,n)] /= L[IDX(j,j,n)];

    // update all remaining rhs
    for (long int k=0; k<Xmat.Cols(); k++)
      for (long int i=j+1; i<n; i++)
	X[IDX(i,k,n)] -= L[IDX(i,j,n)]*X[IDX(j,k,n)];

  }

  // return success
  return 0;
}


// performs forwards substitution on the linear system L*x = b, returning x as a new Mat; 
// leaves L and b untouched
Mat FwdSub(Mat &L, Mat &B) {

  // check that matrix sizes match
  if (L.Rows() != B.Rows() || L.Rows() != L.Cols()) {
    cerr << "FwdSub error, illegal matrix/vector dimensions\n";
    cerr << "  Mat is " << L.Rows() << " x " << L.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols() << endl;
    Mat *X = new Mat(0,0);
    return *X;
  }
  
  // create new Mat for output
  Mat *X = new Mat(L.Rows(),B.Cols());

  // call existing FwdSub routine for computations
  if (FwdSub(L, *X, B) != 0)
    cerr << "FwdSub Warning: error in FwdSub call\n";

  // return result
  return *X;
}


// solves a linear system A*X = B, filling in the input Mat X
int Solve(Mat &Amat, Mat &Xmat, Mat &Bmat) {

  // check that matrix sizes match
  if (Amat.Rows() != Bmat.Rows() || Amat.Rows() != Amat.Cols() ||
      Bmat.Cols() != Xmat.Cols() || Xmat.Rows() != Amat.Rows()) {
    cerr << "Solve error, illegal matrix/vector dimensions\n";
    cerr << "  Mat is " << Amat.Rows() << " x " << Amat.Cols() 
	 << ",  rhs is " << Bmat.Rows() << " x " << Bmat.Cols()
	 << ",  solution is " << Xmat.Rows() << " x " << Xmat.Cols() << endl;
    return 1;
  }
  
  // create temporary variables
  long int i, j, k, p, n;
  double m, tmp, Amax;

  // access the data arrays
  double *A = Amat.get_data();
  double *B = Bmat.get_data();

  // determine maximum absolute entry in A (for singularity check later)
  Amax = Amat.MaxNorm();

  // perform Gaussian elimination to convert A,B to an upper-triangular system
  n = Amat.Rows();
  for (k=0; k<n-1; k++) {   // loop over diagonals

    // find the pivot row p
    p=k;
    for (i=k; i<n; i++)  
      if (fabs(A[IDX(i,k,n)]) > fabs(A[IDX(p,k,n)]))  
	p=i;

    // swap rows in A
    for (j=k; j<n; j++) {
      tmp = A[IDX(p,j,n)];
      A[IDX(p,j,n)] = A[IDX(k,j,n)];
      A[IDX(k,j,n)] = tmp;
    }

    // swap rows in B
    for (j=0; j<Bmat.Cols(); j++) {
      tmp = B[IDX(p,j,n)];
      B[IDX(p,j,n)] = B[IDX(k,j,n)];
      B[IDX(k,j,n)] = tmp;
    }

    // check for nonzero matrix diagonal
    if (fabs(A[IDX(k,k,n)]) < STOL*Amax) {
      cerr << "Solve error: numerically singular matrix!\n";
      return 1;
    }

    // perform elimination using row k
    for (i=k+1; i<n; i++)      // store multipliers in column below pivot
      A[IDX(i,k,n)] /= A[IDX(k,k,n)];
    for (j=k+1; j<n; j++)      // loop over columns of A, to right of pivot 
      for (i=k+1; i<n; i++)    // update rows in column
	A[IDX(i,j,n)] -= A[IDX(i,k,n)]*A[IDX(k,j,n)];
    for (j=0; j<Bmat.Cols(); j++)
      for (i=k+1; i<n; i++)      // update entries in B
	B[IDX(i,j,n)] -= A[IDX(i,k,n)]*B[IDX(k,j,n)];
  }

  // check for singularity at end (only need to check final diagonal entry)
  if (fabs(A[IDX(n-1,n-1,n)]) < STOL*Amax) {
    cerr << "Solve error: numerically singular matrix!\n";
    return 1;
  }

  // check for singularity at end (only need to check final diagonal entry)
  if (fabs(A[IDX(n-1,n-1,n)]) < STOL*Amax) {
    cerr << "Solve error: numerically singular matrix!\n";
    return 1;
  }

  // perform Backwards Substitution on result
  if (BackSub(Amat, Xmat, Bmat) != 0) {
    cerr << "Solve: error in BackSub call\n";
    return 1;
  }

  // return success
  return 0;
}


// solves a linear system A*X = B, returning X as a new Mat
Mat Solve(Mat &A, Mat &B) {

  // check that matrix sizes match
  if (A.Rows() != B.Rows() || A.Rows() != A.Cols()) {
    cerr << "Solve error, illegal matrix/vector dimensions\n";
    cerr << "  Mat is " << A.Rows() << " x " << A.Cols() 
	 << ",  rhs is " << B.Rows() << " x " << B.Cols() << endl;
    Mat *X = new Mat(0,0);
    return *X;
  }

  // create new Mat for output
  Mat *X = new Mat(A.Rows(),B.Cols());

  // call existing Solve routine for computations
  if (Solve(A, *X, B) != 0)
    cerr << "Solve: error in in-place Solve call\n";

  // return result
  return *X;
}


// end of file
