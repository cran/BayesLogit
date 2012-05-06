// -*- mode: c++; fill-column: 80; -*-

//////////////////////////////////////////////////////////////////////
// Jesse Windle - jesse.windle@gmail.com - November, 2011
// See appendix at end of document for additional information.
//////////////////////////////////////////////////////////////////////

#ifndef __MATRIX_FRAME__
#define __MATRIX_FRAME__

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <stdio.h>

using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::ostream;
using std::istream;
using std::stringstream;
using std::ostream;

#ifndef uint
typedef unsigned int uint;
#endif

//////////////////////////////////////////////////////////////////////
                          // MatrixFrame //
//////////////////////////////////////////////////////////////////////

inline bool idxcheck(bool b){
  if (!b) throw std::runtime_error("Index out of bounds.\n");
  return b;
}

inline bool sizecheck(bool b){
  if (!b) throw std::runtime_error("Incompatible dimension.\n");
  return b;
}

inline bool memcheck(bool b){
  if (!b) throw std::runtime_error("Memory overlap.\n");
  return b;
}

class MatrixFrame
{
 protected:
  // I do not want to use *const p since I want to be able to resize my matrices.
  double *p;  // Pointer to the array of matrices.
  uint   nr;  // Number of rows.
  uint   nc;  // Number of columns.
  uint   nm;  // Number of matrices.

  // Matrix Check
  bool allow(uint l) const
    { return (p!=NULL && l < (nr*nc)); }
  bool allow(int  l) const
    { return (p!=NULL && l>=0 && l < (int)(nr*nc)); }
  bool allow(uint r, uint c) const
    { return (p!=NULL && r < nr && c < nc); }
  bool allow(int  r,  int c) const
    { return (p!=NULL && r>=0 && c>=0 && r < (int)nr && c < (int)nc); }

  // Array of Matrix Check
  bool indexok(uint i) const
  { return i < nm; }
  bool indexok( int i) const
  { return (0 <= i && i < (int)nm); }

 public:
  // Constructors, etc.  Do not make a constructor which copies.
  MatrixFrame()
    { p = NULL; nr = 0; nc = 0; nm = 0; }
  MatrixFrame(double *ptr, uint r=1, uint c=1, uint m=1)
    { p = ptr; nr = r; nc = c; nm = m; }
  MatrixFrame(double *ptr,  int r=1,  int c=1,  int m=1)
    { p = ptr; nr = (uint)r; nc = (uint)c; nm = (uint)m; }
  MatrixFrame(const MatrixFrame & M)
    { p = M.p; nr = M.nr; nc = M.nc; nm = M.nm; }

  ~MatrixFrame()
    { p = NULL; }

  // Test equality and assign equality.
  MatrixFrame& operator= (const MatrixFrame &M); // Makes copy.
  bool         operator==(const MatrixFrame &M) const;
  bool         sameframe (const MatrixFrame &M) const;

  // rows(), cols(), mats, area(), vol();
  uint rows() const { return nr; }
  uint cols() const { return nc; }
  uint mats() const { return nm; }
  uint area() const { return nr * nc; }
  uint vol()  const { return nr * nc * nm; }
  int  size() const { return nr * nc * nm; }

  // Matrix Access

  // Returns the (r,c) element of matrix.
  double& operator()(uint r, uint c)
    { idxcheck(allow(r, c)); return p[c * nr + r]; }
  const double& operator()(uint r, uint c) const
    { idxcheck(allow(r, c)); return p[c * nr + r]; }

  // Array of Matrix Access

  // Returns the lth element of array of matrices.
  // I debate whether this is confusing notation.
  const double& operator()(uint l) const
  { idxcheck(l < nr*nc*nm);    return p[l]; }
  double& operator()(uint l)
  { idxcheck(l < nr*nc*nm);    return p[l]; }

  // Returns the (r,c) element of matrix[t].
  double& operator()(uint r, uint c, uint t)
  { idxcheck(indexok(t) && allow(r,c)); return p[t * nc + c * nr + r]; }
  const double& operator()(uint r, uint c, uint t) const
  { idxcheck(indexok(t) && allow(r,c)); return p[t * nc + c * nr + r]; }
  double& get(uint r, uint c=0, uint t=0)
  { idxcheck(indexok(t) && allow(r,c)); return p[t * nc + c * nr + r]; }
  const double& get(uint r, uint c=0, uint t=0) const
  { idxcheck(indexok(t) && allow(r,c)); return p[t * nc + c * nr + r]; }

  // Returns the ith element of the array p.
  double& vec(uint i)
  { idxcheck(i < nr*nc*nm); return p[i]; }
  const double& vec(uint i) const
  { idxcheck(i < nr*nc*nm); return p[i]; }

  // Returns a MatrixFrame pointing to the ith matrix or,
  // if there is one matrix, to the ith column.
  MatrixFrame operator[](uint i)
  { idxcheck(indexok(i)); return MatrixFrame(&p[0+i*area()], nr, nc); }

  // Get the pointer.  Be wary.
  // const double* const getp()
  // { return p; }

  // Array of Matrix Functions.

  void copy(const MatrixFrame& M);     // Copy values.
  // void thincopy(MatrixFrame& M);    // Copy pointer and dimensions.
  MatrixFrame fill(double);            // Fill with value.
  MatrixFrame col(uint c, uint num=1); // The c-th to c+num-1th col.
  MatrixFrame dim(uint r, uint c, uint m=1); // Return a MF with different, compatible dim.

  // Read / Write.
  bool write(      ostream&  os, bool header=0, bool binary=0);
  bool write(const string& file, bool header=0, bool binary=0);

  uint  scan(      istream&  is, bool header=0, bool binary=0);
  uint  scan(const string& file, bool header=0, bool binary=0);
  // bool  readstring(const string& s, bool header=0);

  // Matrix Functions.

  // Fill this matrix from M starting at (r,c).
  void copy(const MatrixFrame& M, uint r, uint c);
  // Copy the matrix M along rows rs and columns cs.
  void copy(const MatrixFrame& M, const MatrixFrame& rs, const MatrixFrame& cs);
  void copy(const MatrixFrame& M, const MatrixFrame& rs, uint c);
  void copy(const MatrixFrame& M, uint r, const MatrixFrame& cs);
  void copy_transpose(const MatrixFrame& M);

  // Set the elements in rows rs and columns cs using the matrix M.
  void set(const MatrixFrame& rs, const MatrixFrame& cs, const MatrixFrame& M);
  void set(const MatrixFrame& rs, uint c, const MatrixFrame& M);
  void set(uint r, const MatrixFrame& cs, const MatrixFrame& M);

}; // MatrixFrame

//////////////////////////////////////////////////////////////////////

#ifndef MF
typedef MatrixFrame MF;
#endif

#ifndef Frame
typedef MatrixFrame Frame;
#endif

//////////////////////////////////////////////////////////////////////
		      // Function Definitions //
//////////////////////////////////////////////////////////////////////

// bool overlap(const MF& a, const MF& b);
// bool hconform(const MF& a, const MF& b);
// uint pconform(const MF& c, const MF& a, const MF& b, char transa='N', char transb='N');
// bool dconform(const MF& a, const MF& b);

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                         IMPLEMENTATION                           //
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
                       // Utility Functions //
//////////////////////////////////////////////////////////////////////

void MatrixFrame::copy(const MatrixFrame& M)
{
  if (this==&M) return;
  idxcheck(nr==M.rows() && nc==M.cols() && nm==M.mats());
  for(uint i = 0; i < vol(); i++) p[i] = M.vec(i);
} // copy

// I am now forcing a person to do this only in a constructor.
// void MatrixFrame::thincopy(MatrixFrame& M)
// {
//   if (this==&M) return;
//   p  = &M(0);
//   nr = M.rows();
//   nc = M.cols();
//   nm = M.mats();
// } // thincopy

MatrixFrame MatrixFrame::col(uint c, uint num)
{
  // Check that these are valid parameters.
  idxcheck(allow((uint)0,c));
  idxcheck(allow((uint)0,c+num-1));
  double *ptr = &operator()(0,c);
  return MatrixFrame(ptr, nr, num);
}

MatrixFrame MatrixFrame::fill(double d)
{
  for(uint i = 0; i < vol(); i++) p[i] = d;
  return *this;
}

MatrixFrame MatrixFrame::dim(uint r, uint c, uint m)
{
  sizecheck (r*c*m==nr*nc*nm);
  return MatrixFrame(p, r, c, m);
}

//////////////////////////////////////////////////////////////////////
                  // Assigment and Test equality //
//////////////////////////////////////////////////////////////////////

/*
  The copy constructor and assignment operator work in two different
  ways.  The copy constructor simply copies the pointer p, and
  dimension info nr and nc into a new object.  These quantities are
  all small and hence there isn't much overhead to this operation.

  The assignment operator works differently.  It copies the contents
  of one MatrixFrame into the contents of another MatrixFrame.  I
  don't want to do this accidentally so I've included a warning.  In
  general, you should not use assignement and instead you should use
  the copy _function_.

  My concern is that one could be using the derived class Matrix and
  write something like c[i] = a[i], which seems to have the more
  intuitive notion of copying the conents of matrix a[i] into matrix
  c[i], but which is still ambiguous since one could just want to copy
  an address.

  To recap, the behavior you should expect:
  shallow copy:
    MatrixFrame mf1(mf2);
    MatrixFrame mf1 = mf2;
  hard copy:
    MatrixFrame mf1, mf2;
    mf1 = mf2;

  Here is something I need to think about.
  Matrix a, b;
  a.col(i) = b.col(j);
  or
  Matrix a = b.col(j);

 */

// Assignment.  See discussion above.
MatrixFrame& MatrixFrame::operator= (const MatrixFrame& M)
{
  //cerr << "Warning: be careful with the assignment operator.\n"
  //     << "MatrixFrame::operator= makes a deep copy.\n";
  if (this==&M) return *this;
  copy(M);
  return *this;
} // operator=

// Test equality.
bool MatrixFrame::operator==(const MatrixFrame& M) const
{
  if(sameframe(M)) return true;
  if(vol() != M.vol()) return false;
  for(uint i = 0; i < vol(); i++) if(M(i)!=operator()(i)) return false;
  return true;
} // operator==

bool MatrixFrame::sameframe(const MatrixFrame& M) const
{
  return (p==&M(0) && nr==M.rows() && nc==M.cols() && nm==M.mats());
}

//////////////////////////////////////////////////////////////////////
			   // Comparison //
//////////////////////////////////////////////////////////////////////

MatrixFrame lt(MatrixFrame c, const MatrixFrame& a, const MatrixFrame& b)
{
  sizecheck(c.vol()==a.vol() && a.vol()==b.vol());
  if(a.sameframe(b)) c.fill(1.0);
  for(uint i = 0; i < a.vol(); i++)
    c(i) = (double) (a(i) < b(i));
  return c;
}

MatrixFrame lteq(MatrixFrame c, const MatrixFrame& a, const MatrixFrame& b)
{
  sizecheck(c.vol()==a.vol() && ( a.vol()%b.vol() )==0 );
  if(a.sameframe(b)) c.fill(0.0);
  for(uint i = 0; i < a.vol(); i++)
    c(i) = (double) (a(i) <= b(i));
  return c;
}

MatrixFrame between(MatrixFrame c, const MatrixFrame a, const MatrixFrame lower, const MatrixFrame upper)
{
  sizecheck(c.vol()==a.vol() && a.vol()==lower.vol() && lower.vol()==upper.vol());
  for(uint i = 0; i < a.vol(); i++)
    c(i) = (double) ( lower(i) <= a(i) && a(i) <= upper(i) );
  return c;
}

// MatrixFrame within(MatrixFrame c, const MatrixFrame a, double lower, double upper)
// {
//   sizecheck(c.vol()==a.vol() && a.vol()==lower.vol() && lower.vol()==upper.vol());
//   for(uint i = 0; i < a.vol(); i++)
//     c(i) = (double) ( lower <= a(i) && a(i) <= upper );
//   return c;
// }

//////////////////////////////////////////////////////////////////////
		    // Input / Output Operators //
//////////////////////////////////////////////////////////////////////

// Output matrix to stream.  This is done in a column major way so
// that a human reading the output will see an array of transposed
// matrices.

ostream& operator<<(ostream& os, MatrixFrame M)
{
  M.write(os, false, false);
  return os;
}

// Read in data from a string using scan.

MatrixFrame& operator<<(MatrixFrame& M, const string& s)
{
  stringstream ss(s);
  M.scan(ss, false, false);
  return M;
}

//////////////////////////////////////////////////////////////////////
			  // Read / Write //
//////////////////////////////////////////////////////////////////////

// Writes a matrix.  You may chose to include a header, which are the
// dimensions of the array of matrices.

bool MatrixFrame::write(std::ostream& os, bool header, bool binary)
{
  if (!os) return false;
  // Write binary.
  if(binary){
    if(header){
      os.write((char*) &nr, sizeof(nr));
      os.write((char*) &nc, sizeof(nc));
      os.write((char*) &nm, sizeof(nm));
    }
    for (uint i = 0; i < vol(); i++)
      os.write((char*) &p[i], sizeof(double));
  }
  // Write human.
  if(!binary){
    if(header)
      os << nr << " " << nc << " " << nm << "\n";
    for(uint k = 0; k < nm; k++){
      for(uint j = 0; j < nc; j++){
	for(uint i = 0; i < nr; i++){
	  os << operator()(i,j,k) << " ";
	}
	os << "\n";
      }
      if ((k+1) != nm) os << "\n";
    }
  }
  os.flush();
  return true;
} // write

bool MatrixFrame::write(const string& file, bool header, bool binary)
{
  std::ofstream ofs(file.c_str());
  if (!ofs) return false;
  return write(ofs, header, binary);
} // write

// Reads in data from a string of values until the end of the stream
// or the end of the array of matrices is reached.  You are alerted if
// you do not read in enough data to fill the array.

uint MatrixFrame::scan( std::istream& is, bool header, bool binary)
{
  // Tell us if something is wrong.
  if (!is || is.eof())  return 0;

  uint i = 0; // The nubmer of items read.

  // Read binary.
  if(binary){
    if(header){
      uint r,c,m;
      is.read((char*) &r, sizeof(nr));
      is.read((char*) &c, sizeof(nc));
      is.read((char*) &m, sizeof(nm));
      // sizecheck(vol() == r*c*m); // A somewhat strict condition.
    }
    while(!is.eof() && i < vol())
      is.read((char*) &p[i++], sizeof(double));
  }
  // Write human.
  if(!binary){
    if(header){
      uint r,c,m;
      is >> r >> c >> m;
      // sizecheck(vol() == r*c*m); // A somewhat strict condition.
    }
    while(!is.eof() && i < vol())
      { is >> p[i++]; ws(is); } // ws extracts intermediate white space.
                                // Needed in case the stream is padded by white space.
  }
  // Warnings:
  if (i != vol())
    Rprintf( "In scan: Number of items read (%i) different \
                     than number of elements in Matrix.\n", i);
  if (!is.eof())
    Rprintf( "In scan: Did not reach end of file.\n");

  return i;
} // scan

uint MatrixFrame::scan(const string& file, bool header, bool binary)
{
  std::ifstream ifs(file.c_str());
  if(!ifs){
    Rprintf( "Cannot read file %s.\n", file.c_str());
    return 0;
  }
  return scan(ifs, header, binary);
} // read

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//                                                                  //
//     MATRIX OPERATIONS - IE THINGS THAT ONLY USE FIRST MATRIX     //
//                                                                  //
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
                  // Conformity and Overlap Check //
//////////////////////////////////////////////////////////////////////

// Check if MF a and MF b overlap in memory.

bool overlap(const MF& a, const MF& b)
{
  // Get the low and the high pointer.
  const MatrixFrame low  = &a(0) < &b(0) ? a : b;
  const MatrixFrame high = &a(0) < &b(0) ? b : a;
  // Check if there is overlap.  Do we need to subtract 1?  Yes.
  return (&low(0) + low.area() - 1) < &high(0) ? false : true;
} // overlap


// Hadamard conform.  This is not commutative.  Checks if the size of
// a is equally divided by the size of b.

bool hconform(const MF& a, const MF& b)
{
  return ( a.area() % b.area() )==0;
} // hconform

// Matrix product conform.  Checks if c = op(a) * op(b) is valid.
// Returns 0 if invalid and the matching dimension if valid.
uint pconform(const MF& c, const MF& a, const MF& b, char transa='N', char transb='N')
{
  uint opa_rows = transa=='T' ? a.cols() : a.rows();
  uint opa_cols = transa=='T' ? a.rows() : a.cols();
  uint opb_rows = transb=='T' ? b.cols() : b.rows();
  uint opb_cols = transb=='T' ? b.rows() : b.cols();
  bool conform  = (opa_cols==opb_rows) && (c.rows()==opa_rows) && (c.cols()==opb_cols);
  if (!conform) {
    Rprintf("c_rows: %u\n", c.rows());
    Rprintf("c_cols: %u\n", c.cols());
    Rprintf("opa_rows: %u\n", opa_rows);
    Rprintf("opa_cols: %u\n", opa_cols);
    Rprintf("opb_rows: %u\n", opb_rows);
    Rprintf("opb_cols: %u\n", opb_cols);
  }
  return conform ? opa_cols : 0;
} // pconform

// Matrices are the same dimension.
bool dconform(const MF& a, const MF& b)
{
  return ( a.rows()==b.rows() && a.cols()==b.cols() );
}

//////////////////////////////////////////////////////////////////////
			    // COPYING //
//////////////////////////////////////////////////////////////////////

// Need to check for alisasing in memory!!!

// This matrix is smaller or the same area as M.
// Fill this matrix from M starting at (r,c).
void MatrixFrame::copy(const MatrixFrame& M, uint r, uint c)
{
  sizecheck( (r + nr) <= M.rows() && (c + nc) <= M.cols() );
  for(uint j = 0; j < nc; j++)
    for(uint i = 0; i < nr; i++)
      operator()(i,j) = M(r+i, c+j);
} // copy

void MatrixFrame::copy(const MatrixFrame& M, const MatrixFrame& rs, const MatrixFrame& cs)
{
  sizecheck(rs.area()==rows() && cs.area()==cols());
  // Should check min and max of rs and cs to make sure you are in bounds.
  // I suppose this is checked by indexing.
  for(uint j = 0; j < cs.area(); j++){
    for(uint i = 0; i < rs.area(); i++){
      operator()(i,j) = M(rs(i), cs(j));
    }
  }
} // copy

void MatrixFrame::copy(const MatrixFrame& M, const MatrixFrame& rs, uint c)
{
  sizecheck(rs.area()==rows() && 1==cols());
  for(uint i = 0; i < rs.area(); i++){
    operator()(i,0) = M(rs(i), c);
  }
} // copy

void MatrixFrame::copy(const MatrixFrame& M, uint r, const MatrixFrame& cs)
{
  sizecheck(1==rows() && cs.area()==cols());
  for(uint j = 0; j < cs.area(); j++){
    operator()(0,j) = M(r, cs(j));
  }
} // copy

void MatrixFrame::copy_transpose(const MatrixFrame& M)
{
  sizecheck(nr==M.cols() && nc==M.rows() && nm==M.mats());
  for(uint k = 0; k < nm; ++k){
    for(uint j = 0; j < nc; ++j){
      for(uint i = 0; i < nr; ++i){
	get(i,j,k) = M.get(j,i,k);
      }
    }
  }
}

void MatrixFrame::set(const MatrixFrame& rs, const MatrixFrame& cs, const MatrixFrame& M)
{
  sizecheck(rs.area()==M.rows() && cs.area()==M.cols());
  for(uint j = 0; j < cs.area(); j++){
    for(uint i = 0; i < rs.area(); i++){
      operator()(rs(i),cs(j)) = M(i, j);
    }
  }
}

void MatrixFrame::set(const MatrixFrame& rs, uint c, const MatrixFrame& M)
{
  sizecheck(rs.area()==M.area());
  for(uint i = 0; i < rs.area(); i++){
    operator()(rs(i), c) = M(i);
  }
}

void MatrixFrame::set(uint r, const MatrixFrame& cs, const MatrixFrame& M)
{
  sizecheck(cs.area()==M.area());
  for(uint j = 0; j < cs.area(); j++){
    operator()(r, cs(j)) = M(j);
  }
}

//////////////////////////////////////////////////////////////////////
                      // Hadamard Operations //
//////////////////////////////////////////////////////////////////////

// A Hadamard operation is an operation done element wise.
// Notationally, given op in {*,+,/,-} these functions perform
// something like

//    a[i] op= b[i % b.area()]
//    c[i] = alpha * a[i] op b[i % b.area()] + beta c[i].

// Thus the elements of b are recycled when area(b) < area(a).  We
// require that area(b) be a multiple of area(a).

// The functions available are h<op>eq(a, b, sc) and h<op>(c, a, b,
// alpha=0.0, beta=1.0) where <op> is prod, sum, div, or sub.

// Hadamard Operation Equals (HOPEQ) a op= sc * b
#define HOPEQ(NAME, OPEQ)						\
  MatrixFrame NAME(MF a, const MF& b, double sc=1.0)			\
  {                                                                     \
    memcheck(!overlap(a,b));                                              \
    sizecheck(hconform(a,b));                                              \
    uint barea = b.area();						\
    for(uint i = 0; i < a.area(); i++)					\
      a(i) OPEQ sc * b(i % barea);					\
    return a;								\
  }                                                                     \
  MatrixFrame NAME(MF a, double b)					\
  {                                                                     \
    for(uint i = 0; i < a.area(); i++)					\
      a(i) OPEQ b;							\
    return a;								\
  }                                                                     \
  MatrixFrame& operator OPEQ(MF& a, const MF& b)			\
  {                                                                     \
    memcheck(!overlap(a,b));                                              \
    sizecheck(hconform(a,b));                                              \
    uint barea = b.area();						\
    for(uint i = 0; i < a.area(); i++)					\
      a(i) OPEQ b(i % barea);						\
    return a;								\
  }									\

HOPEQ(hprodeq, *=) HOPEQ(hsumeq, +=)
HOPEQ(hdiveq,  /=) HOPEQ(hsubeq, -=)

#undef HOPEQ

// Hadamard Operation (HOP) c = a op sc * b
#define HOP(NAME, OP)                                                   \
  void NAME(MF c, const MF& a, const MF& b, double alpha=1.0, double beta=0.0)	\
  {                                                                     \
    bool okay = (!overlap(c,b) && !overlap(a,b)) ||			\
      (!overlap(c,a) && !overlap(c,b));					\
    memcheck(okay);							\
    sizecheck(hconform(a,b));                                              \
    sizecheck(c.area()==a.area());						\
    uint barea = b.area();						\
    for(uint i = 0; i < c.area(); i++)					\
      c(i) = alpha * a(i) OP b(i % barea) + beta * c(i);		\
  }                                                                     \
  void NAME(MF c, const MF& a, double b)                                      \
  {                                                                     \
    sizecheck(c.area()==a.area());						\
    for(uint i = 0; i < c.area(); i++)					\
      c(i) = a(i) OP b;                                                 \
  }                                                                     \
  void NAME(MF c, double b, const MF& a)					\
  {                                                                     \
    sizecheck(c.area()==a.area());						\
    for(uint i = 0; i < c.area(); i++)					\
      c(i) = b OP a(i);                                                 \
  }                                                                     \

HOP(hprod, *) HOP(hsum, +)
HOP(hdiv,  /) HOP(hsub, -)

#undef HOP

//////////////////////////////////////////////////////////////////////
                      // Operation along rows //
//////////////////////////////////////////////////////////////////////

// Sometime we want to take an operation "along rows" i.e. the
// operation a(i,j) op= b(j % b.area()) where op may be *,+,/,-.  The
// functions to do this are <op>onrow where <op> is prod, sum, div, or
// sub.  We require that a.cols() to be a multiple of area(b).

#define ROWOP(NAME, OPEQ)			 \
  MatrixFrame NAME(MF a, MF b)			 \
  {						 \
    memcheck(!overlap(a, b));			 \
    sizecheck(a.cols()%b.area()==0);		 \
    uint barea = b.area();			 \
    for(uint j = 0; j < a.cols(); j++)		 \
      for(uint i = 0; i < a.rows(); i++)	 \
        a(i,j) OPEQ b(j % barea);		 \
    return a;					 \
  }						 \

ROWOP(prodonrow, *=) ROWOP(sumonrow, +=)
ROWOP(divonrow,  /=) ROWOP(subonrow, -=)

#undef ROWOP

//////////////////////////////////////////////////////////////////////
			// Statistics //
//////////////////////////////////////////////////////////////////////

double sum(const MF& a)
{
  double total = 0.0;
  for(uint i = 0; i < a.vol(); i++)
    total += a(i);
  return total;
}

double mean(const MF& a)
{
  return sum(a) / a.vol();
}

//////////////////////////////////////////////////////////////////////
			 // BLAS / LAPACK //
//////////////////////////////////////////////////////////////////////

/*
  See BLAS / LAPACK documentation at netlib.org.  The Fortran source
  code is very regular.  The perl script BLAStoC.perl will take a BLAS
  subroutine/function file and convert it to the necessary C code.
 */

//////////////////////////////////////////////////////////////////////
		       // WRAPPER TO FORTRAN //
//////////////////////////////////////////////////////////////////////

extern "C" {

  // BLAS LEVEL 1 //

  void daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);

  double ddot_(int* N, double* DX, int* INCX, double* DY, int* INCY);

  // BLAS LEVEL 3 //

  void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA, double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);

  void dtrmm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* ALPHA, double* A, int* LDA, double* B, int* LDB);

  void dtrsm_(char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* ALPHA, double* A, int* LDA, double* B, int* LDB);

  // LAPACK //

  void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);

  void dposv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB, int* INFO);

  void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);

}

//////////////////////////////////////////////////////////////////////
		    // MatrixFrame BLAS WRAPPER //
//////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------//
// y = alpha x + y.

void daxpy(int n, double da, double* dx, int incx, double* dy, int incy)
{ daxpy_(&n, &da, dx, &incx, dy, &incy); }

void axpy(double alpha, MF x, MF y)
{
  sizecheck(x.rows()==y.rows() && x.cols()==1 && y.cols()==1);
  daxpy((int)x.rows(), alpha, &x(0), 1, &y(0), 1);
}

//------------------------------------------------------------------//
// x'y

double ddot(int n, double* dx, int incx, double* dy, int incy)
{ return ddot_(&n, dx, &incx, dy, &incy); }

double dot(MF x, MF y)
{
  sizecheck(x.rows()==y.rows() && x.cols()==1 && y.cols()==1);
  return ddot(x.rows(), &x(0), 1, &y(0), 1);
}

//------------------------------------------------------------------//
// c = alpha op(a) * op(b) + beta c.

void dgemm(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b, int ldb, double beta, double* c, int ldc)
{ dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

void gemm(MF c, MF a, MF b, char ta='N', char tb='N', double alpha=1.0, double beta=0.0)
{
  memcheck(!overlap(c,a) && !overlap(c,b));
  // Get the dimensionality information we need.
  int cnr = (int)c.rows(); int cnc = (int)c.cols();
  int anr = (int)a.rows(); int bnr = (int)b.rows();
  int k   = (int)pconform(c, a, b, ta, tb);
  // Make sure things conform.
  sizecheck(k!=0);
  dgemm(ta, tb, cnr, cnc, k, alpha, &a(0), anr, &b(0), bnr, beta, &c(0), cnr);
} // gemm

//------------------------------------------------------------------//
// b = alpha op(a) * b  OR  b = alpha b * op(a) where a is triangular.

void dtrmm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
{ dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

void trmm(MF a, MF b, char uplo, char side='L', char ta='N', char diag='N', double alpha=1.0)
{
  memcheck(!overlap(a,b));
  // This checks that a is square and that the product conforms.
  uint k = side=='L' ? pconform(b, a, b, ta, 'N') : pconform(b, b, a, 'N', ta);
  sizecheck(k!=0);
  dtrmm(side, uplo, ta, diag, b.rows(), b.cols(), alpha, &a(0), a.rows(), &b(0), b.rows());
} // trmm

//------------------------------------------------------------------//
// Solve x:  op(a) x = alpha b  OR  x op(a) = alpha b, a triangular.
// i.e: x = alpha inv(op(a)) b  OR  x = alpha b inv(op(a)).
// The solution is overwriten into b.

void dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, double* a, int lda, double* b, int ldb)
{ dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

void trsm(MF a, MF b, char uplo, char side='L', char ta='N', char diag='N', double alpha=1.0)
{
  memcheck(!overlap(a,b));
  // This checks that a is square and that the product conforms.
  uint k = side=='R' ? pconform(b, a, b, ta, 'N') : pconform(b, b, a, 'N', ta);
  sizecheck(k!=0);
  dtrsm(side, uplo, ta, diag, b.rows(), b.cols(), alpha, &a(0), a.rows(), &b(0), b.rows());
} // trsm

//////////////////////////////////////////////////////////////////////
		  // MATRIX FRAME LAPACK WRAPPER //
//////////////////////////////////////////////////////////////////////

// Solve a general linear system, ax = b for x.

void dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int& info)
{ dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); }

int gesv(MF a, MF b)
{
  memcheck(!overlap(a, b));       // No overlap in memory.
  sizecheck(pconform(b, a, b)!=0); // a is square and b conforms.
  int info;
  std::vector<int> ipiv(a.rows());
  dgesv(a.rows(), b.cols(), &a(0), a.rows(), &ipiv[0], &b(0), b.rows(), info);
  return info;
}

// Shorthand.
int solve(MF a, MF b)
{
  return gesv(a, b);
}

//------------------------------------------------------------------//
// Solves ax = b for x where a is sym. pos. def.  Note: the lower (or
// upper) portion of A is overwritten with the Cholesky decomposition.

void dposv(char uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info)
{ dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info); }

int posv(MF a, MF b, char uplo)
{
  memcheck(!overlap(a,b));
  sizecheck(pconform(b, a, b)!=0);
  int info;
  dposv(uplo, a.rows(), b.cols(), &a(0), a.rows(), &b(0), b.rows(), info);

  if (info != 0) {
    Rprintf("Error in posv: info = %i\n", info);
    throw std::runtime_error("aborted in posv\n");
  }

  return info;
}

//------------------------------------------------------------------//
// Cholesky Decomposition

void dpotrf(char uplo, int n, double* a, int lda, int& info)
{ dpotrf_(&uplo, &n, a, &lda, &info); }

int potrf(MF a, char uplo)
{
  sizecheck(a.rows()==a.cols());
  int info = 0;
  dpotrf(uplo, a.rows(), &a(0), a.rows(), info);
  return info;
}

int chol(MF a, char uplo='L')
{
  return potrf(a, uplo);
}

//------------------------------------------------------------------//


//////////////////////////////////////////////////////////////////////
			  // END OF CODE //
//////////////////////////////////////////////////////////////////////

#endif // MATRIX_FRAME

//////////////////////////////////////////////////////////////////////
			    // APPENDIX //
//////////////////////////////////////////////////////////////////////

/*********************************************************************

  The goal of the matrixMatrixFrame class and the Matrix class is to
  represent arrays of matrices.  An array with only one matrix can be
  thought of simply as a matrix.  We had three goals in mind when
  creating this clss: 1) keep it simple/transparent, 2) make it easy
  to use in MCMC, and 3) make it easy to talk to R.

  Regarding (1), there is a tradeoff between ease of understanding the
  code and ease of calculations.  Eigen is a great Matrix package, but
  it is hard to understand what exactly is going on.  I don't
  understand expression templates, but they make calculations nice and
  easy.  Since this will be used in MCMC simulations we wanted to know
  exactly what is going on under the hood, which means that this code
  is more understandable but that you will pay a cost when expressing
  your computations.

  Regarding (2), we hope that an array of Matrices will be sufficient
  for most MCMC algorithms.  It is possible that one may need an array
  of an array of matrices, such as an Matrix Normal DLM within a Gibbs
  sampler, but we haven't ecountered that too often in our MCMC
  experience.

  Regarding (3), you should be able to pass the pointer to some data
  from R directly to the MatrixFrame class to set up your data
  structure.

  You can think of the MatrixFrame class as simply a wrapper to
  BLAS/LAPACK, and you can think of the Matrix class as a simple data
  container.  That's all there is too it.  In that sense, MatrixFrame
  is an "interface".  The matrix operations are ALWAYS performed on
  the first matrix in the array.  If you want to perform an operation
  on different matrix in the array use the [] operator.

  Below you will find that the class is split into two categories:
  those parts of the class that have to do with the array of matrices
  and those parts of the class that have to do with the first matrix
  in the array, i.e. the matrix on which we are operating.

  In general, we need to be careful with this class.  There are lots
  of opportunities to break things since the class possess a pointer
  to some data, which the class itself did not create.

  Everything is COLUMN MAJOR for compatibility with Fortran.

  Idea: You could template this to the dimension of the array.

*********************************************************************/
