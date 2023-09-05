#ifndef INCLUDED_VECTOR_H
#define INCLUDED_VECTOR_H               //!< Prevent multiple inclusion
#include<vector>
#include<assert.h>
#include<config.h>
/*!\file
  \brief  Declaration/definition of Vector and Matrix classes.

  Also declarations of functions that perform operations on complex vectors.

*/
#include <string.h>
#include <comms/sysfunc_cps.h>
#include <util/data_types.h>
#include <util/vector_asm.h>
#include <util/verbose.h>
#include <util/omp_wrapper.h>

#define VEC_INLINE

#if 0
#define NOINLINE_MACRO __attribute((noinline))
#else
#define NOINLINE_MACRO 
#endif

CPS_START_NAMESPACE

const int OMP_CUTOFF=100;
class Vector; 	// forward declaration
class Matrix;

//------------------------------------------------------------------
// Declarations of some genaral c-style functions that perform
// operations on vectors. For these functions there exists 
// optimized assembly code.
//------------------------------------------------------------------
extern "C" 
{
    //! vector copy; b = a
//void moveMem(void *b, const void *a, int len); 
inline void moveMem(void *b, const void *a, int len) 
{
#undef PROFILE
#ifdef PROFILE
    double time  = -dclock();
#endif
    memcpy(b, a, len); 
#ifdef PROFILE
    time += dclock();
    print_flops("","moveMem",len,time);
#endif
}

//void moveFloat(Float *b, const Float *a, int len); 
inline void moveFloat(Float *b, const Float *a, int len) {

    if(len> OMP_CUTOFF){
#pragma omp parallel for
    for(int i =0;i<len;i++) b[i] = a[i];
    } else {
        memcpy(b, a, len*sizeof(Float)); 
    }
}

    //! 3x3 complex matrix multiplication; c = ab 
#ifndef VEC_INLINE
void mDotMEqual(IFloat* c, const IFloat* a, const IFloat* b);
#else
inline void mDotMEqual(IFloat* c, const IFloat* a, const IFloat* b)
{
    *c      = *a      * *b      - *(a+1)  * *(b+1)    +
    	      *(a+2)  * *(b+6)  - *(a+3)  * *(b+7)    +
    	      *(a+4)  * *(b+12) - *(a+5)  * *(b+13);
    *(c+1)  = *a      * *(b+1)  + *(a+1)  * *b        +
    	      *(a+2)  * *(b+7)  + *(a+3)  * *(b+6)    +
    	      *(a+4)  * *(b+13) + *(a+5)  * *(b+12);

    *(c+2)  = *a      * *(b+2)  - *(a+1)  * *(b+3)    +
    	      *(a+2)  * *(b+8)  - *(a+3)  * *(b+9)    +
    	      *(a+4)  * *(b+14) - *(a+5)  * *(b+15);
    *(c+3)  = *a      * *(b+3)  + *(a+1)  * *(b+2)    +
    	      *(a+2)  * *(b+9)  + *(a+3)  * *(b+8)    +
    	      *(a+4)  * *(b+15) + *(a+5)  * *(b+14);

    *(c+4)  = *a      * *(b+4)  - *(a+1)  * *(b+5)    +
    	      *(a+2)  * *(b+10) - *(a+3)  * *(b+11)   +
    	      *(a+4)  * *(b+16) - *(a+5)  * *(b+17);
    *(c+5)  = *a      * *(b+5)  + *(a+1)  * *(b+4)    +
    	      *(a+2)  * *(b+11) + *(a+3)  * *(b+10)   +
    	      *(a+4)  * *(b+17) + *(a+5)  * *(b+16);

    *(c+6)  = *(a+6)  * *b      - *(a+7)  * *(b+1)    +
    	      *(a+8)  * *(b+6)  - *(a+9)  * *(b+7)    +
    	      *(a+10) * *(b+12) - *(a+11) * *(b+13);
    *(c+7)  = *(a+6)  * *(b+1)  + *(a+7)  * *b        +
    	      *(a+8)  * *(b+7)  + *(a+9)  * *(b+6)    +
    	      *(a+10) * *(b+13) + *(a+11) * *(b+12);

    *(c+8)  = *(a+6)  * *(b+2)  - *(a+7)  * *(b+3)    +
    	      *(a+8)  * *(b+8)  - *(a+9)  * *(b+9)    +
    	      *(a+10) * *(b+14) - *(a+11) * *(b+15);
    *(c+9)  = *(a+6)  * *(b+3)  + *(a+7)  * *(b+2)    +
    	      *(a+8)  * *(b+9)  + *(a+9)  * *(b+8)    +
    	      *(a+10) * *(b+15) + *(a+11) * *(b+14);

    *(c+10) = *(a+6)  * *(b+4)  - *(a+7)  * *(b+5)    +
    	      *(a+8)  * *(b+10) - *(a+9)  * *(b+11)   +
    	      *(a+10) * *(b+16) - *(a+11) * *(b+17);
    *(c+11) = *(a+6)  * *(b+5)  + *(a+7)  * *(b+4)    +
    	      *(a+8)  * *(b+11) + *(a+9)  * *(b+10)   +
    	      *(a+10) * *(b+17) + *(a+11) * *(b+16);

    *(c+12) = *(a+12) * *b      - *(a+13) * *(b+1)    +
    	      *(a+14) * *(b+6)  - *(a+15) * *(b+7)    +
    	      *(a+16) * *(b+12) - *(a+17) * *(b+13);
    *(c+13) = *(a+12) * *(b+1)  + *(a+13) * *b        +
    	      *(a+14) * *(b+7)  + *(a+15) * *(b+6)    +
    	      *(a+16) * *(b+13) + *(a+17) * *(b+12);

    *(c+14) = *(a+12) * *(b+2)  - *(a+13) * *(b+3)    +
    	      *(a+14) * *(b+8)  - *(a+15) * *(b+9)    +
    	      *(a+16) * *(b+14) - *(a+17) * *(b+15);
    *(c+15) = *(a+12) * *(b+3)  + *(a+13) * *(b+2)    +
    	      *(a+14) * *(b+9)  + *(a+15) * *(b+8)    +
    	      *(a+16) * *(b+15) + *(a+17) * *(b+14);

    *(c+16) = *(a+12) * *(b+4)  - *(a+13) * *(b+5)    +
    	      *(a+14) * *(b+10) - *(a+15) * *(b+11)   +
    	      *(a+16) * *(b+16) - *(a+17) * *(b+17);
    *(c+17) = *(a+12) * *(b+5)  + *(a+13) * *(b+4)    +
    	      *(a+14) * *(b+11) + *(a+15) * *(b+10)   +
    	      *(a+16) * *(b+17) + *(a+17) * *(b+16);
}
#endif

    //! CK: 3x3 complex matrix multiplication with complex conjugate on first matrix; c = Conj(a)b 
void mStarDotMEqual(IFloat* c, const IFloat* a, const IFloat* b);

  //! CK: 3x3 complex matrix multiplication with complex conjugate on second matrix; c = a Conj(b) 
void mDotMStarEqual(IFloat* c, const IFloat* a, const IFloat* b);

  //! CK: 3x3 complex matrix multiplication with complex conjugate on second matrix; c = Conj(a) Conj(b) 
void mStarDotMStarEqual(IFloat* c, const IFloat* a, const IFloat* b);

    //! 3x3 complex matrix multiplication and sum; c += ab
void mDotMPlus(IFloat* c, const IFloat* a, const IFloat* b); 

 //! CK: 3x3 complex matrix multiplication and sum; c += Conj(a)b
void mStarDotMPlus(IFloat* c, const IFloat* a, const IFloat* b); 

 //! CK: 3x3 complex matrix multiplication and sum; c += a Conj(b)
void mDotMStarPlus(IFloat* c, const IFloat* a, const IFloat* b); 

//! CK: 3x3 complex matrix multiplication and sum; c += Conj(a) Conj(b)
void mStarDotMStarPlus(IFloat* c, const IFloat* a, const IFloat* b); 

    //! 3x3 complex matrix times vector; y = Mx
void uDotXEqual(IFloat* y, const IFloat* m, const IFloat* x); 

    //! vector scalar product; a.b
IFloat dotProduct(const IFloat *a, const IFloat *b, int);

    //! vector addition; a += b
#ifndef VEC_INLINE
void vecAddEquVec(IFloat *a, const IFloat *b, int); 	
#else 
inline void vecAddEquVec(IFloat *a, const IFloat *b, int len)
{
#pragma omp parallel for
    for(int i = 0; i < len; ++i) {
    	a[i] += b[i];
    }
}
#endif

    //! vector subtraction; a -= b
void vecMinusEquVec(IFloat *a, const IFloat *b, int);  

inline void vecMinusEquVecSingle(IFloat *a, const IFloat *b, int len)
{
    for(int i = 0; i < len; ++i) {
        *a++ -= *b++;
    }
}

    //! vector negation; a = -b
void vecNegative(IFloat *a, const IFloat *b, int); 	

    //! set all elements to zero
void vecZero(IFloat *a, int size);

/*!
  \param a The vector to be multiplied
  \param b The real scalar
  \param len The length of the vectors.
  \post \a a is the multiplied vector.
 */
inline void vecTimesEquFloat(IFloat *a, IFloat b, size_t len) {
#pragma omp parallel for
    for(size_t i = 0; i < len; ++i) {
    	a[i] *= b;
    }
}

inline void vecAddEquFloat(IFloat *a, IFloat b, size_t len) {
#pragma omp parallel for
    for(size_t i = 0; i < len; ++i) {
        a[i] += b;
    }
}

inline void vecTimesEquFloatSingle(IFloat *a, IFloat b, size_t len)
{
    for(size_t i = 0; i < len; ++i) {
    	*(a+i) *= b;
    }
}

void vecTimesComplex(IFloat *a,
                     IFloat re,
                     IFloat im,
                     const IFloat *c,
                     int len);

    //! real scalar times vector multiplication; a = c*b
void vecEqualsVecTimesEquFloat(IFloat *a, IFloat *b, IFloat c, int); // 


    //! vector linear combination; a = bc+d
void fTimesV1PlusV2(IFloat *a, IFloat b, const IFloat *c,
			const IFloat *d, int size); 	

inline void fTimesV1PlusV2Single(IFloat *a, IFloat b, const IFloat *c,
        const IFloat *d, int len)
{
    for(int i = 0; i < len; ++i) {
        *a++ = b * *c++ + *d++;
    }
}

    //! vector linear combination; a = bc-d
void fTimesV1MinusV2(IFloat *a, IFloat b, const IFloat *c,
			const IFloat *d, int size);    

    //! complex vector scalar product; a.b
void compDotProduct(IFloat *c_r, IFloat *c_i, 
        	    const IFloat *a, const IFloat *b, int);

    //! complex vector linear combination; a = bc+d
void cTimesV1PlusV2(IFloat *a, IFloat b_re, IFloat b_im, const IFloat *c,
                    const IFloat *d, int size);      

    //! Not implemented on qcdsp
void cTimesV1MinusV2(IFloat *a, IFloat b_re, IFloat b_im, const IFloat *c,
	             const IFloat *d, int size);      // A = b*C-D

    //! matrix linear combination; a = 1-bc
void oneMinusfTimesMatrix(IFloat *a, IFloat b, const IFloat *c, int n);     

}


//------------------------------------------------------------------
// Declarations of some genaral c-style functions that perform
// operations on vectors. For these functions there is
// no optimized assembly.
//------------------------------------------------------------------

//! Multiplication of complex vector by matrix and addition; y += Mx
void uDotXPlus(IFloat* y, const IFloat* u, const IFloat* x);

//! Multiplication of complex vector by matrix and subtraction; y -= Mx
void uDotXMinus(IFloat* y, const IFloat* u, const IFloat* x);

//! Multiplication of complex vector by hermitian conjugate matrix and summation; y += M^dagger x
void uDagDotXEqual(IFloat* y, const IFloat* u, const IFloat* x);

//! Multiplication of complex vector by hermitian conjugate matrix; y = M^dagger x
void uDagDotXPlus(IFloat* y, const IFloat* u, const IFloat* x);

//------------------------------------------------------------------
// Declarations of some genaral c-style functions that compute
// re/im parts of the su3 characters of various representations
// of su3 matrices.
//------------------------------------------------------------------

extern IFloat reChar6(IFloat *p) ;
extern IFloat imChar6(IFloat *p) ;
extern IFloat reChar8(IFloat *p) ;
extern IFloat reChar10(IFloat *p) ;
extern IFloat imChar10(IFloat *p) ;


//------------------------------------------------------------------
//! The rank of the matrices represented by the Matrix class
//------------------------------------------------------------------
enum{COLORS=3}; 

//------------------------------------------------------------------
//! A class of general 3x3 complex matrices.
//------------------------------------------------------------------
class Matrix
{
    Float u[2*COLORS*COLORS];	// The matrix

    friend class Vector;
    static IFloat inv3;

  public:
    // CREATORS
    //! General constructor; no initialisation.
    Matrix() {}

    //! Initialisation to real multiple of the unit matrix.
    //------------------------------------------------------------------
    /*!
      The diagonal matrix elements (0,0), (1,1) and (2,2) are set to the real
      number \a c; All other elements are zero.  
      \param c The diagonal matrix element
    */
    Matrix(IFloat c) {*this = c;}

    //! Initialisation to complex multiple of the unit matrix.
    //------------------------------------------------------------------
    /*!
      The diagonal matrix elements (0,0), (1,1) and (2,2) are set to the complex
      number \a c; All other elements are zero.  
      \param c The diagonal matrix element
    */
    Matrix(const Complex& c) {*this = c;}

    //! Copy constructor
    //------------------------------------------------------------------
    /*!
      The matrix is initialised as a copy of the matrix \a m.
      \param m The initialising matrix.
    */
    Matrix(const Matrix& m) {
        memcpy(u, m.u, sizeof(Float) * COLORS * COLORS * 2);
    }

    //! Assignment to real multiple of the unit matrix.
    //------------------------------------------------------------------
    /*!
      The diagonal matrix elements (0,0), (1,1) and (2,2) are set to the real
      number \a c; All other elements are zero.  
      \param c The diagonal matrix element
    */
    Matrix& operator=(IFloat c) {
        this->ZeroMatrix();
        u[0] = u[8] = u[16] = c;
        return *this;
    }

    //! Assignment to complex multiple of the unit matrix.
    //------------------------------------------------------------------
    /*!
      The diagonal matrix elements (0,0), (1,1) and (2,2) are set to the complex
      number \a c; All other elements are zero.  
      \param c The diagonal matrix element
    */
    Matrix& operator=(const Complex& c) {
        this->ZeroMatrix();
        u[0] = u[8] = u[16] = c.real();
        u[1] = u[9] = u[17] = c.imag();
        return *this;
    }
    //! Overloaded assignment
    /*! \a m should not alias this matrix */
    Matrix& operator=(const Matrix& m) {
        if(this != &m) 
            memcpy(u, m.u, sizeof(Float) * COLORS * COLORS * 2);
        return *this;
    }

    // MANIPULATORS
    //! Adds a matrix \a m to this matrix.
    /*!
      \param m The matrix to be added.
      \return The matrix sum.
    */
    Matrix& operator+=(const Matrix& m) {
        for(int i = 0; i < COLORS * COLORS * 2; ++i)
            u[i] += m.u[i];
        return *this;
    }

    Matrix& operator*=(const Matrix& m) {
      Matrix tmp(*this);
      mDotMEqual((IFloat *)u, (IFloat *) tmp.u, (IFloat *) m.u);
      //this->DotMEqual(tmp,m);
      return *this;
    }


    //! Adds a real scalar multiple of the unit matrix to this one.
    /*!
      \param c The real scalar multiple
      \return The matrix sum
    */
    Matrix& operator+=(IFloat c)
     { u[0] += c;  u[8] += c;  u[16] += c;  return *this; }

    //! Subtracts a matrix \a m to this matrix.
    /*!
      \param m The matrix to be subtracted.
      \return The matrix difference.
    */
    Matrix& operator-=(const Matrix& m)
    { vecMinusEquVecSingle((IFloat *)u, (IFloat *)m.u, COLORS*COLORS*2);
      return *this; }

    //! Subtracts a real scalar multiple of the unit matrix from this one.
    /*!
      \param c The real scalar multiple
      \return The matrix difference
    */
    Matrix& operator-=(IFloat c)
      { u[0] -= c;  u[8] -= c;  u[16] -= c;  return *this; }

    //! Multiplies this matrix by a real scalar.
    /*!
      \param c The real scalar
      \return The multiplied matrix
    */
    Matrix& operator*=(IFloat c) {
        for(int i = 0; i < COLORS * COLORS * 2; ++i)
            u[i] *= c;
        return *this;
    }
    //! Multiplies this matrix by a complex scalar.
    /*!
      \param c The complex scalar
      \return The multiplied matrix
    */

#if 1
    Matrix& operator*=(const Complex &c) {
      for(int i = 0; i < COLORS * COLORS; ++i){
	int reidx = 2*i; int imidx = reidx+1;
	Float rev = u[reidx];
	u[reidx] = u[reidx]*c.real() - u[imidx]*c.imag();
	u[imidx] = rev*c.imag() + u[imidx]*c.real();
      }
      return *this;
    }
#else
    // Added by Hantao
    Matrix &operator*=(const Complex &c) {
        Complex *uc = (Complex *)u;
        for(int i = 0; i < COLORS * COLORS; ++i) {
            uc[i] *= c;
        }
        return *this;
    }
#endif

    Matrix operator+(const Matrix &m)const {
        Matrix tmp(*this);
        tmp += m;
        return tmp;
    }

    Matrix operator-(const Matrix &m)const {
        Matrix tmp(*this);
        tmp -= m;
        return tmp;
    }

    Matrix operator*(const Matrix &m)const {
        Matrix tmp;
        tmp.DotMEqual(*this, m);
        return tmp;
    }

    Matrix operator*(const Complex &c) {
      Matrix tmp(*this);
      for(int i = 0; i < COLORS * COLORS; ++i){
	int reidx = 2*i; int imidx = reidx+1;
	Float rev = tmp.u[reidx];
	tmp.u[reidx] = u[reidx]*c.real() - u[imidx]*c.imag();
	tmp.u[imidx] = rev*c.imag() + u[imidx]*c.real();
      }
      return tmp;
    }


     //! Assignment to matrix product; \a ab
     /*!
         \param a the matrix \a a
	 \param b the matrix \a b
	 \post This matrix is the matrix product \a ab
      */
    void DotMEqual(const Matrix& a, const Matrix& b)
    { mDotMEqual((IFloat *)u, (IFloat *) a.u, (IFloat *) b.u);}

     //! Assignment to matrix product; \a ab
     /*!
         \param a the matrix \a a
	 \param b the matrix \a b
	 \post The matrix product \a ab  is added to  this matrix.
      */
    void DotMPlus(const Matrix& a, const Matrix& b)
    { mDotMPlus((IFloat *)u, (IFloat *)a.u, (IFloat *)b.u);}
       //u += a.u * b.u

    //! Assignment to Matrix transpose.
    void Trans(const IFloat* m);

    //! Assignment to matrix transpose.
    /*!
      \param m A matrix.
      \post This matrix is the transpose of \a m.
      
      \a m must not be an alias of this matrix/
    */
    void Trans(const Matrix& m)
        { Trans((const IFloat *)(m.u)); }

    //! Assignment to Matrix complex conjugate.
    void Conj(const IFloat* m);

    //! Assignment to matrix complex conjugate.
    /*!
      \param m A matrix.
      \post This matrix is the complex conjugate of \a m.
      
      \a m must not be an alias of this matrix/
    */
    void Conj(const Matrix& m)
        { Conj((const IFloat *)(m.u)); }

    //! Hermitian conjugate.
    void Dagger(const Matrix& m)
    	{ Dagger((const IFloat *)&m); }
    void Dagger(){
       Matrix dag;
       dag.Dagger(*this);
       *this = dag;
    }
    void Transpose(const IFloat* m);
    void Transpose();

    //! Determine matrix trace
    Complex Trace() const{
      return Complex(u[0]+u[8]+u[16],u[1]+u[9]+u[17]);
    }
      
    //! Assignment to hermitian conjugate.
    /*!
      \param m A matrix.
      \post This matrix is the hermitian conjugate of \a m.

      \a a must not be an alias of this matrix
    */
#ifndef VEC_INLINE
    void Dagger(const IFloat* m);

    //! Not what you might think.
    void TrLessAntiHermMatrix(const Matrix& this_dag);
//    void TrLessAntiHermMatrix();
    void TrLessAntiHermMatrix(){
       Matrix dag;
       dag.Dagger(*this);
       this->TrLessAntiHermMatrix(dag);
    }
#else 
//#define TAH_INLINE
/*!
  \param a A linear array representation of a 3x3 complex matrix, such that 
  real part of the (i,j) element is at array position [6i+2j] 
  and the imaginary part of the (i,j) element is at array position [6i+2j+1].
  \post This matrix is the hermitian conjugate of \a m.

  \a a must not be an alias of this matrix.
 */
inline void Dagger(const IFloat* a)
{
    u[0]  = a[0];   u[1]  = -a[1];
    u[6]  = a[2];   u[7]  = -a[3];
    u[12] = a[4];   u[13] = -a[5];
    u[2]  = a[6];   u[3]  = -a[7];
    u[8]  = a[8];   u[9]  = -a[9];
    u[14] = a[10];  u[15] = -a[11];
    u[4]  = a[12];  u[5]  = -a[13];
    u[10] = a[14];  u[11] = -a[15];
    u[16] = a[16];  u[17] = -a[17];
}

/*!
  \param dag A matrix \a A.
  \post This matrix is set to\n
  <em>1/2(M-A) - 1/6 Trace M-A)</em>
  \n where \a M is the original value of this matrix.
*/
inline void TrLessAntiHermMatrix(const Matrix& dag)
{
    // get 1/2(A - dag(A)) =  1/2A - dag(1/2A)
    *this -= dag;

    IFloat *p = (IFloat *)u;
    vecTimesEquFloatSingle(p, 0.5, 18);

    IFloat c = inv3 * (*(p+1) + *(p+9) + *(p+17));
    *(p+1) -= c;
    *(p+9) -= c;
    *(p+17) -= c;
}

inline void TrLessAntiHermMatrix()
{

    IFloat *p = (IFloat *)u;
    *p = *(p+8) = *(p+16)=0.;
    IFloat tmp = 0.5*(p[2] - p[6]);
    p[2]=tmp; p[6] = -tmp;
    tmp = 0.5*(p[3] + p[7]);
    p[3]=tmp; p[7] = tmp;
    tmp = 0.5*(p[4] - p[12]);
    p[4]=tmp; p[12] = -tmp;
    tmp = 0.5*(p[5] + p[13]);
    p[5]=tmp; p[13] = tmp;
    tmp = 0.5*(p[10] - p[14]);
    p[10]=tmp; p[14] = -tmp;
    tmp = 0.5*(p[11] + p[15]);
    p[11]=tmp; p[15] = tmp;

    IFloat c = inv3 * (*(p+1) + *(p+9) + *(p+17));
    p[1] -= c;
    p[9] -= c;
    p[17] -= c;
}
#endif

    //! Assignment to tensor product of vectors.
    void Cross2(const Vector& v1, const Vector& v2);

    //! Assignment to an traceless antihermitian matrix.
    void AntiHermMatrix(const IFloat *a);
    	// a points to an array of 8 real numbers
	    // *this = i \lamda^i * a^i
	    // \lambda^i are 8 Gellmann matrices

    //! Force this matrix to be an SU(3) matrix.
    void Unitarize(void);

    //! Only do the last step in Unitarize.
    void Construct3rdRow(void);

    //! Project this matrix onto SU(3) according to its polar decomposition
    //! Added by Thomas Dumitrescu 06/2006
    int ProjSU3(void);

    //! Assignment to a unit matrix.
    //------------------------------------------------------------------
    /*!
      \post This matrix is a 3x3 unit matrix.
    */
    void UnitMatrix(void) {
        this->ZeroMatrix();
        u[0] = u[8] = u[16] = 1.;
    }

    //! Assignment to a zero matrix.
    //------------------------------------------------------------------
    /*!
      \post This matrix is a 3x3 zero matrix.
    */
    void ZeroMatrix(void) {
        for(int i = 0; i < 18; ++i) {
            u[i] = 0.;
        }
    }

    //! Assignment to a negated matrix.
    /*!
      \param m A matrix.
      \post This matrix has the value \a -m.
    */
    void NegMatrix(const Matrix& m) {
        for(int i = 0; i < COLORS * COLORS * 2; ++i)
            u[i] = -m.u[i];
    }

    //! Assignment to the matrix linear combination 1-xm
    /*!
      \param x A real scalar factor
      \param m A matrix
      \post This matrix has the value 1-xm.
    */
    void OneMinusfTimesM(IFloat x, const Matrix& m)
    	{ oneMinusfTimesMatrix((IFloat *)u, x, (IFloat *)&m,
	  COLORS*COLORS*2); }


#ifndef VEC_INLINE
    // ACCESSORS
    //! Write access.
    Complex& operator()(int i, int j);
    //! Read access.
    const Complex& operator()(int i, int j) const;
#else
Complex& operator()(int i, int j)
{ return ((Complex*)u)[i*COLORS+j]; }
const Complex& operator()(int i, int j) const
{ return ((Complex*)u)[i*COLORS+j]; }
#endif

    //! Write access.
    /*!
      \param i A number between 0 and 8
      \return The ([i - i mod 3]/3, i mod 3) matrix element

      Should this method not be private?
    */
    inline Complex& operator[](int i) { return ((Complex*)u)[i]; }
    //! Read access.
    /*!
      \param i A number between 0 and 8
      \return The ([i - i mod 3]/3, i mod 3) matrix element
    */
    inline const Complex& operator[](int i) const { return ((Complex*)u)[i]; }
    inline IFloat elem(int i) { return u[i]; }
    //! Read access.
    /*!
      \param i A number between 0 and 17
      \return  element of the array
    */
    //! The determinant.
    void Det(IFloat* c) const;

    //! Returns the real part of the trace.
    IFloat ReTr() const {
        return u[0] + u[8] + u[16];
    }

    //! Returns the trace.
    Complex Tr() const {
        return ((Complex*)u)[0] + ((Complex*)u)[4] + ((Complex*)u)[8];
    }

    //! -1/2 times the trace of the square.
    IFloat NegHalfTrSquare() const;

    //! The deviation of this matrix from unitarity
    IFloat ErrorSU3() const;

    /*! Returns the SU(3) matrix norm, defined by
      ||X||^2 = -2 trace X^2
    */
    // !!FIXME: Why does it calculate this?
    IFloat norm() const {
        Matrix x2;
        x2.DotMEqual(*this, *this);
        return -2.0*x2.ReTr();
    }
    IFloat norm2() const {
        IFloat *m = (IFloat*)&u[0];
        return dotProduct(m, m, 18);
    }
    // SU(3) Characters

    Complex Char3() const { return Tr() ; } ;

    Complex Char6() const ;

    Complex Char8() const ;

    Complex Char10() const ;

    void FTimesV1PlusV2(Float fb, Matrix *c, Matrix *d, int len)
	{ fTimesV1PlusV2((IFloat *)&u, IFloat(fb), (IFloat *)c, 
	                         (IFloat *)d, len*18); }
    Float Norm()
    {
      Float sum=0.;
      for(int i=0; i<2*COLORS*COLORS; i++) sum +=u[i]*u[i];
      return sum;
    }


};


inline static Matrix Transpose(const Matrix &m){
  Matrix out;
  out.Trans(m);
  return out;
}

//Added by CK
inline static Complex Trace(const Matrix &a, const Matrix &b){
  //Mapping is i*3 + j  
  Complex out(0.0);
  //a(0,0)*b(0,0) + a(0,1)*b(1,0) + a(0,2)*b(2,0)
  out += a[0]*b[0] + a[1]*b[3] + a[2]*b[6]; 
  //a(1,0)*b(0,1) + a(1,1)*b(1,1) + a(1,2)*b(2,1)
  out += a[3]*b[1] + a[4]*b[4] + a[5]*b[7];
  //a(2,0)*b(0,2) + a(2,1)*b(1,2) + a(2,2)*b(2,2)
  out += a[6]*b[2] + a[7]*b[5] + a[8]*b[8];
  return out;
}




//------------------------------------------------------------------
//! A class implementing a general 3 component complex vector.
/*!
  This is a schizophrenic class. It is really designed to be a class
  of complex 3-vectors, and many methods carry out operations on just
  such an object; \e e.g. the overloaded binary operators, the
  matrix-vector multiplications and the normalisation and orthogonalisation
  methods. However, some methods, those with take an argument \c int \a len,
  are really wrappers for functions operating on 1-dimensional floating
  point arrays of any length. They are meant to be used with an array of
  Vectors: the first Vector in the array operates not only on its own data
  but on that of all the other objects in the array by assuming that it is
  at the beginning of a contiguous floating point array. For the sake of
  sanity the argument \a len should be a multiple of 6.
*/
//------------------------------------------------------------------
class Vector
{
    Float v[2*COLORS];	// Vector

    friend class Matrix;

  public:
    // CREATORS
    Vector() {}

    //! Overloaded assignment
    /*! \a x should not alias this matrix */
    Vector& operator=(const Vector& x)
#if 1
    { for(int i=0;i<COLORS*2;i++) v[i] = x.v[i];
#else
    { moveMem(v, x.v, COLORS*2*sizeof(Float)); 
#endif
      return *this; }

    // MANIPULATORS
    //! Multiplies this vector by a real scalar.
    /*!
      \param p The real scalar
      \return The multiplied vector
    */
    Vector& operator*=(IFloat p)
    { vecTimesEquFloatSingle((IFloat *)v, p, COLORS*2); 
    return *this; }

    //! Adds a vector \a m to this vector.
    /*!
      \param m The vector to be added.
      \return The vector sum.
    */
    Vector& operator+=(const Vector& x)
    { vecAddEquVec((IFloat *)v, (IFloat *)x.v, COLORS*2);
      return *this; }

    //! Subtracts a vector \a m to this vector.
    /*!
      \param m The vector to be subtracted.
      \return The vector difference.
    */
    Vector& operator-=(const Vector& x)
    { vecMinusEquVecSingle((IFloat *)v, (IFloat *)x.v, COLORS*2);
      return *this; }

    //! Assignment to matrix-vector product.
    /*!
      \param m A matrix.
      \param x a vector
      \post This vector is takes the value Mx
    */
    void DotXEqual(const Matrix& m, const Vector& x)
    { uDotXEqual((IFloat *)v, (IFloat *) m.u, (IFloat *) x.v); }
       // v = m.u * x.v, m should be in CRAM, x MUST be in DRAM */

    //! Normalisation
    void Normalize(void);

    //! Orthogonalisation
    void Orthogonalize(const Vector& p1);

    //! Zeroing a color vector 
    /*!
      added by Sam 1/9/2006 to implement disconnected F.T
    */
    void Zero()
    {
      for(int i=0; i<2*COLORS; i++) v[i]=0;
    }

    //! simple element access (as a Complex)
    /*!
      added by Sam 1/9/2006 to implement disconnected F.T
    */
    Complex& operator[](int i) 
    {
      return *((Complex*)(v+2*i));
    }

    //--------------------------------------------------------------
    // Functions that act on arrays of vectors of general length.
    // The array of vectors is treated as an array of IFloating 
    // numbers pointed to by &v and having length len.
    // This set of functions does not really fit the way
    // Vector is currently defined (as an array with re/im and 
    // color indeces only. It extends the notion of Vector to
    // a general array of IFloating numbers.
    //--------------------------------------------------------------

    //! Assignment to another vector
    /*!
      \param b Another vector
      \param len The number of real numbers in the vectors.
      \post This vector  = \a b

      \a b should not alias this vector.
    */
    void CopyVec(const Vector *b, int len)
#if 1
    { moveFloat((Float *)&v, (const Float *)b, len); }
#else
    { moveMem(&v, b, len*sizeof(Float)); }
#endif

    //! Square norm.
    /*!
      \param len The number of real numbers in the vectors.
      \return The square norm of this vector.
    */
    Float NormSqNode(size_t len)
    {return Float( dotProduct((IFloat *)&v, (IFloat *)&v, len) ); }

    //! Square norm with global sum.
    Float NormSqGlbSum(size_t len);
    Float NormSqGlbSum4D(size_t len);
    //! Print the vector content to the screen
    void Print(int len);
    
    //! The real part of the dot product
    /*!
      \param b Another vector
      \param len The number of real numbers in the vectors.
      \return The real part of the dot product (v,b).
    */
    Float ReDotProductNode(const Vector *b, int len)
    {return Float( dotProduct((IFloat *)&v, (IFloat *)b, len) ); }

    //! The real part of the dot product with global sum.
    Float ReDotProductGlbSum(const Vector *b, int len);
    Float ReDotProductGlbSum4D(const Vector *b, int len);

    void NormSqArraySliceSum(Float *f_out, const int size, const int dir);
    //!< Sum the square norms of vectors in 3-dim slices.

    void SliceArraySum(Float *sum, const Float *f_in, const int dir);
    //!< Sum an array of Floats on a 4-dim lattice in 3-dim slices.

    void SliceArraySumFive(Float *sum, const Float *f_in, const int dir);
    //!< Sum an array of Floats on a 5-dim lattice in 4-dim slices.

    //! Assign vector to zero.
    /*!
      \param len The number of real numbers in the vectors.
      \post This vector has the value 0.
    */
    void VecZero(int len)
      {vecZero((IFloat*)&v, len);}

    //! Assignment to a negated vector.
    /*!
      \param b A vector.
      \param len The number of real numbers in the vectors.
      \post This vector has the value \a -b.
    */
    void VecNegative(const Vector *b, int len)
	{vecNegative((IFloat *)&v, (IFloat *)b, len);}

    //! Multiplication by a real scalar
    /*!
      \param fb The real scalar
      \param len The number of real numbers in the vectors.
      \post This vector is multiplied by \a fb
    */
    void VecTimesEquFloat(const Float &fb, size_t len)
    {vecTimesEquFloat((IFloat *)&v, IFloat(fb), len);}

    //! Multiplication by a real scalar
    /*!
      \param u The input vector
      \param fb The real scalar
      \param len The number of real numbers in the vectors.
      \post This vector is multiplied by \a fb
    */
    void VecEqualsVecTimesEquFloat(const Vector *u, const Float &fb, int len)
    {vecEqualsVecTimesEquFloat((IFloat *)&v, (IFloat*)u, IFloat(fb), len);}

    //! Addition of another vector
    /*!
      \param b Another vector
      \param len The number of real numbers in the vectors.
      \post \a b is added to this vector.
    */
    void VecAddEquVec(const Vector *b, int len)
    { vecAddEquVec((IFloat *)&v, (IFloat *)b, len);}

    //! Subtraction of another vector
    /*!
      \param b Another vector
      \param len The number of real numbers in the vectors.
      \post \a b is subtracted from this vector.
    */
    void VecMinusEquVec(const Vector *b, int len)
    { vecMinusEquVec((IFloat *)&v, (IFloat *)b, len);}

    //! Assignment of the linear combination  fb * c + d
    /*!
      \param fb A real scalar
      \param c A vector
      \param d A vector      
      \param len The number of real numbers in the vectors.
      \post \a  This vector takes the value fb * c + d
    */
//    void FTimesV1PlusV2(const Float &fb, const Vector *c,
//			const Vector *d, int len)
     void FTimesV1PlusV2(Float fb, Vector *c, Vector *d, int len)
    { fTimesV1PlusV2((IFloat *)&v, IFloat(fb), (IFloat *)c, 
                        (IFloat *)d, len); }
     void FTimesPlusVec(Float fb, Vector *c, int len){
		IFloat *a_f = (IFloat *)&v;
		IFloat *c_f = (IFloat *)&(c->v);

#pragma omp parallel for default(shared)
    for(int i = 0; i < len; ++i) {
    	*(a_f+i) += fb* *(c_f+i);
    }

	}


    //! Assignment of the linear combination  fb * c - d
    /*!
      \param fb A real scalar
      \param c A vector
      \param d A vector      
      \param len The number of real numbers in the vectors.
      \post  This vector takes the value fb * c - d
    */
    void FTimesV1MinusV2(const Float &fb, const Vector *c,
			 const Vector *d, int len)
    { fTimesV1MinusV2((IFloat *)&v, IFloat(fb), (IFloat *)c, 
                        (IFloat *)d, len); }
    //! The dot product with another vector
    /*!
      \param b Another vector
      \param len The number of real numbers in the vectors.
      \return The dot product of this vector with b.
    */
    Complex CompDotProductNode(const Vector *b, int len)
    {
       IFloat c_r, c_i;
       compDotProduct(&c_r, &c_i, (IFloat *)&v, (IFloat *)b, len);
       return Complex(c_r,c_i);
     }

    //! The dot product with another vector, with global sum
     Complex CompDotProductGlbSum(const Vector *b, int len);
     Complex CompDotProductGlbSum4D(const Vector *b, int len);
  
    //! Assignment of the linear combination  fb * c + d
    /*!
      \param fb A complex scalar
      \param c A vector
      \param d A vector      
      \param len The number of real numbers in the vectors.
      \post  This vector takes the value fb * c + d
    */
    void CTimesV1PlusV2(const Complex &fb, const Vector *c,
                         const Vector *d, int len)
    { cTimesV1PlusV2((IFloat *)&v, real(fb), imag(fb), (IFloat *)c, 
	              (IFloat *)d, len); }
 
    //! Assignment of the linear combination  fb * c - d
    /*!
      \param fb A complex scalar
      \param c A vector
      \param d A vector      
      \param len The number of real numbers in the vectors.
      \post \a  This vector takes the value fb * c - d
    */
    void CTimesV1MinusV2(const Complex &fb, const Vector *c,
                           const Vector *d, int len)
    { cTimesV1MinusV2((IFloat *)&v, real(fb), imag(fb), (IFloat *)c, 
	                (IFloat *)d, len); }

    void print(const char *name,size_t f_size){
//       Float *v_p = (Float*)v;
       Float sum=NormSqGlbSum(f_size);
       VRB.Result("",name,"%0.12g %0.12g %0.12g %0.12g %0.12g %0.12g norm=%0.12g\n",
		v[0],v[1],v[2],v[3],v[4],v[5],sum);
    }


};

  inline void vaxpy3(Vector *res,Float *scale,Vector *mult,Vector
*add, int ncvec){
  fTimesV1PlusV2((IFloat *)res, (IFloat)*scale, (IFloat *)mult,
    (IFloat *)add, ncvec*6);
}
  inline void vaxpy3_m(Matrix *res,Float *scale,Matrix *mult,Matrix
*add, int ncvec){
  fTimesV1PlusV2((IFloat *)res, (IFloat)*scale, (IFloat *)mult,
    (IFloat *)add, ncvec*6);
}

#if 0
inline void moveFloattofloat NOINLINE_MACRO (float *out, Float * in, size_t f_size)
{
  Float  sum=0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t i = 0; i < f_size; i++) {
      out[i]=in[i];
//    flt = (float) in[i];
//    out[i] = flt;
      sum +=out[i]*out[i];
  }
  glb_sum(&sum);
  VRB.Result("","moveFloattofloat()","norm=%e\n",sum);
};

inline void movefloattoFloat NOINLINE_MACRO (Float * out, float *in, size_t f_size)
{
//  float flt;
  Float  sum=0.;
#pragma omp parallel for reduction(+:sum)
  for (size_t i = 0; i < f_size; i++) {
      out[i]=in[i];
//    flt = in[i];
//    out[i] = (Float) flt;
      sum +=out[i]*out[i];
  }
  glb_sum(&sum);
  VRB.Result("","moveFloattofloat()","norm=%e\n",sum);
};
#endif

template < typename AFloat, typename BFloat >
  void compDotProduct (std::vector < Float > &result,
                       const std::vector < AFloat * >a,
                       const std::vector < BFloat * >b, size_t len)
{

  const char *fname="compDotProduct()";
  size_t a_size = a.size ();
  size_t b_size = b.size ();
  size_t c_size = len;

  int a_step = 8;
  int b_step = 8;
  int c_step = 8;
  if (a_size < a_step) a_step = a_size;
  if (b_size < b_step) b_step = b_size;
  if (c_size < c_step ) c_step = c_size;
  VRB.Debug("",fname,"sizes= %d %d %d step = %d %d %d result.size()=%d\n", 
	a_size,b_size,c_size,a_step,b_step,c_step,result.size());fflush(stdout);
  int nthr=0;
#pragma omp parallel 
{
   nthr = omp_get_num_threads();
}
//  assert ((len % c_step) == 0);
//  assert ((a_size % a_step) == 0);
//  assert ((b_size % b_step) == 0);
  result.resize (2 * a_size * b_size, 0.);
  Float *result_p = result.data();
//  exit(-42);
//#pragma omp parallel for reduction(+:re,im)
  for (size_t i = 0; i < a_size; i += a_step)
    for (size_t j = 0; j < b_size; j += b_step)
#pragma omp parallel for
      for (size_t k = 0; k < c_size; k += c_step) {
        int ii_end = a_step; if((a_size-i)<a_step) ii_end = a_size-i;
        int jj_end = b_step; if((b_size-j)<b_step) jj_end = b_size-j;
        int kk_end = c_step; if((c_size-k)<c_step) kk_end = c_size-k;
//  printf("%s:end= %d %d %d \n", fname, ii_end,jj_end,kk_end);fflush(stdout);
//OMP4( parallel for )
        for (size_t kk = 0; kk < kk_end; kk++) {
          std::vector < Float > re (a_step * b_step*2, 0);
          size_t ind_k = 2 * (k + kk);
          for (size_t ii = 0; (ii < ii_end); ii++) {
            AFloat *a_p = a[i + ii];
            for (size_t jj = 0; (jj < jj_end); jj++) {
              BFloat *b_p = b[j + jj];
//              printf("%s: %d %d: (%e %e) (%e %e)\n", fname, 
//              if(!ind_k) VRB.Result("",fname,"%d %d: (%e %e) (%e %e)\n",
//		i+ii,j+jj, *a_p,*(a_p+1),*b_p, *(b_p+1));fflush(stdout);
              re[2*(ii + a_step * jj)] +=
                *(a_p + ind_k) * *(b_p + ind_k) + *(a_p + ind_k + 1) * *(b_p + ind_k + 1);
              re[1+2*(ii + a_step * jj)] +=
                *(a_p + ind_k) * *(b_p + ind_k + 1) - *(a_p + ind_k + 1) * *(b_p + ind_k);
            }
          }
#if 1
          for (size_t ii = 0; ii < ii_end; ii++)
            for (size_t jj = 0; jj < jj_end ; jj++){
//OMP4( critical )
#pragma omp critical
            {
              size_t ind_ij = (i + ii) + a_size * (j + jj);
              result_p[2 * ind_ij] += re[2*(ii + a_step * jj)];
              result_p[2 * ind_ij + 1] += re[1+2*(ii + a_step * jj)];
            }
          }
#endif
        }
      }

}

CPS_END_NAMESPACE
#endif
