/* physics system */

#include <alg/wilson_matrix.h>

CPS_START_NAMESPACE

#ifndef INLINE_WILSON_MATRIX

#ifndef USE_C11
inline void cmad( Rcomplex& x, const Rcomplex& y, const Rcomplex& z )
{
#if 1
  x += y*z;
#else
  x.real()+=y.real()*z.real();
  x.real()-=y.imag()*z.imag();
  x.imag()+=y.imag()*z.real();
  x.imag()+=y.real()*z.imag();
#endif
}

inline void cmeq( Rcomplex& x, const Rcomplex& y, const Rcomplex& z )
{
#if 1
  x = y*z;
#else
  x.real() =y.real()*z.real();
  x.real()-=y.imag()*z.imag();
  x.imag() =y.imag()*z.real();
  x.imag()+=y.real()*z.imag();
#endif
}
#else
inline void cmad( Rcomplex& x, const Rcomplex& y, const Rcomplex& z )
{ x += y * z; }
inline void cmeq( Rcomplex& x, const Rcomplex& y, const Rcomplex& z )
{ x = y * z; }
#endif

WilsonMatrix& eq_mult( WilsonMatrix& xmat,
		       const WilsonMatrix& amat,
		       const WilsonMatrix& bmat )
{
  const Rcomplex* a(amat.ptr());
  const Rcomplex* b(bmat.ptr());
  Rcomplex* xoff(xmat.ptr());
  Rcomplex const *point;
  for (int i1=0;i1<12;++i1)
    {
      point = b;
      const Rcomplex& aval(*a);
      cmeq(xoff[0] ,aval, point[0]);
      cmeq(xoff[1] ,aval, point[1]);
      cmeq(xoff[2] ,aval, point[2]);
      cmeq(xoff[3] ,aval, point[3]);
      cmeq(xoff[4] ,aval, point[4]);
      cmeq(xoff[5] ,aval, point[5]);
      cmeq(xoff[6] ,aval, point[6]);
      cmeq(xoff[7] ,aval, point[7]);
      cmeq(xoff[8] ,aval, point[8]);
      cmeq(xoff[9] ,aval, point[9]);
      cmeq(xoff[10],aval, point[10]);
      cmeq(xoff[11],aval, point[11]);
      a++;
      point+=12;
      for (int i3=1;i3<12;++i3)
	{
	  const Rcomplex& aval(*a);
	  cmad(xoff[0] ,aval, point[0]);
	  cmad(xoff[1] ,aval, point[1]);
	  cmad(xoff[2] ,aval, point[2]);
	  cmad(xoff[3] ,aval, point[3]);
	  cmad(xoff[4] ,aval, point[4]);
	  cmad(xoff[5] ,aval, point[5]);
	  cmad(xoff[6] ,aval, point[6]);
	  cmad(xoff[7] ,aval, point[7]);
	  cmad(xoff[8] ,aval, point[8]);
	  cmad(xoff[9] ,aval, point[9]);
	  cmad(xoff[10],aval, point[10]);
	  cmad(xoff[11],aval, point[11]);
 	  a++;
	  point+=12;
	}
      xoff+=12;
    }
  return xmat;
}
#endif

CPS_END_NAMESPACE
