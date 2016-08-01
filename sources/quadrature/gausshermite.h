#ifndef "GAUSSHERMITE_H"
#define "GAUSSHERMITE_H"

# include <cstdlib>
# include <cstdio>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

void gen_hermite_compute ( int order, double alpha, double x[], double w[] );
void gen_hermite_handle ( int order, double alpha, int option, string output );
void gen_laguerre_compute ( int order, double alpha, double xtab[], 
  double weight[] );
void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] );
void gen_laguerre_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_gamma ( double x );
double r8_huge ( );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

void gen_hermite_compute ( int order, double alpha, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule.
//
//  Discussion:
//
//    The integral to be approximated has the form:
//
//      Integral ( -oo < x < +oo ) x^ALPHA exp(-x^2) f(x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double ALPHA, the parameter.
//
//    Output, double X[ORDER], W[ORDER], the abscissas and weights
//    for the requested generalized Gauss-Hermite rule.
//
{
  double alpha_laguerre;
  double arg;
  int i;
  int order_laguerre;
  double *w_laguerre;
  double *x_laguerre;

  if ( order == 1 )
  {
    arg = ( alpha + 1.0 ) / 2.0;
    x[0] = 0.0;
    w[0] = r8_gamma ( arg );
    return;
  }

  if ( ( order % 2 ) == 0 ) 
  {
    order_laguerre = order / 2;
    alpha_laguerre = ( alpha - 1.0 ) / 2.0;
  }
  else
  {
    order_laguerre = ( order - 1 ) / 2;
    alpha_laguerre = ( alpha + 1.0 ) / 2.0;
  }
  
  w_laguerre = new double[order_laguerre];
  x_laguerre = new double[order_laguerre];

  gen_laguerre_compute ( order_laguerre, alpha_laguerre, 
    x_laguerre, w_laguerre );

  if ( ( order % 2 ) == 0 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i];
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+i] = 0.5 * w_laguerre[i];
    }
  }
  else if ( ( order % 2 ) == 1 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    x[order_laguerre] = 0.0;
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+1+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i] / x_laguerre[order_laguerre-1-i];
    }

    arg = ( alpha + 1.0 ) / 2.0;
    w[order_laguerre] = r8_gamma ( arg );
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre] = w[order_laguerre] - w_laguerre[i] / x_laguerre[i];
    }

    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+1+i] = 0.5 * w_laguerre[i] / x_laguerre[i];
    }
  }
  delete [] w_laguerre;
  delete [] x_laguerre;

  return;
}
//****************************************************************************80

void gen_hermite_handle ( int order, double alpha, int option, string output )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_HERMITE_HANDLE computes a generalized Gauss-Hermite rule and outputs it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double ALPHA, the parameter.
//
//    Input, int OPTION.
//    * 0, the integral has the form 
//      Integral ( -oo < x < +oo ) x^ALPHA exp(-x^2) f(x) dx
//    * 1, the integral has the form 
//      Integral ( -oo < x < +oo )                   f(x) dx.
//
//    Input, string OUTPUT_FILE, specifies the output.
//    * "C++'", print as C++ code.
//    * "F77", print as FORTRAN77 code.
//    * "F90", print as FORTRAN90 code.
//    * "MAT", print as MATLAB code.
//    * file,  write files "file_w.txt", "file_x.txt", "file_r.txt" 
//      defining weights, abscissas, and region.
// 
{
  int i;
  string output_r;
  string output_w;
  string output_x;
  double r[2];
  double *w;
  double *x;

  r[0] = - r8_huge ( );
  r[1] =   r8_huge ( );

  x = new double[order];
  w = new double[order];

  gen_hermite_compute ( order, alpha, x, w );
//
//  Modify weights if requested.
//
  if ( option != 0 )
  {
    if ( ( order % 2 ) == 0 )
    {
      for ( i = 0; i < order; i++ )
      {
        w[i] = exp ( x[i] * x[i] ) * w[i] / pow ( r8_abs ( x[i] ), alpha );
      }
    }
    else
    {
      for ( i = 0; i < ( order / 2 ); i++ )
      {
        w[i] = exp ( x[i] * x[i] ) * w[i] / pow ( r8_abs ( x[i] ), alpha );
      }
      for ( i = ( order / 2 ) + 1; i < 2 * ( order / 2 ) + 1; i++ )
      {
        w[i] = exp ( x[i] * x[i] ) * w[i] / pow ( r8_abs ( x[i] ), alpha );
      }
    }
  }
//
//  Print out the rule.
//
  if ( output == "C++" )
  {
    cout << "//\n";
    cout << "//  Weights W, abscissas X and range R\n";
    cout << "//  for a generalized Gauss-Hermite quadrature rule\n";
    cout << "//  ORDER = " << order << "\n";
    cout << "//  ALPHA = " << alpha << "\n";
    cout << "//\n";

    if ( option == 0 )
    {
      cout << "//  OPTION = 0, Standard rule:\n";
	  cout << "//    Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx\n";
	  cout << "//    is to be approximated by\n";
	  cout << "//    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    else
    {
      cout << "//  OPTION = 1, modified rule:\n";
      cout << "//    Integral ( -oo < x < +oo ) f(x) dx\n";
      cout << "//    is to be approximated by\n";
      cout << "//    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    cout << "//\n";

    for ( i = 0; i < order; i++ )
    {
      cout << "  w[" << i << "] = " 
           << setprecision(16) << w[i] << ";\n";
    }
    cout << "\n";
    for ( i = 0; i < order; i++ )
    {
      cout << "  x[" << i << "] = " 
           << setprecision(16) << x[i] << ";\n";
    }
    cout << "\n";
    for ( i = 0; i < 2; i++ )
    {
      cout << "  r[" << i << "] = " << r[i] << ";\n";
    }
   }
   else if ( output == "F77" )
   {
    cout << "c\n";
    cout << "c  Weights W, abscissas X and range R\n";
    cout << "c  for a generalized Gauss-Hermite quadrature rule\n";
    cout << "c  ORDER = " << order << "\n";
    cout << "c  ALPHA = " << alpha << "\n";
    cout << "c\n";
    if ( option == 0 )
    {
      cout << "c  OPTION = 0, Standard rule:\n";
	  cout << "c    Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx\n";
	  cout << "c    is to be approximated by\n";
	  cout << "c    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    else
    {
      cout << "c  OPTION = 1, modified rule:\n";
      cout << "c    Integral ( -oo < x < +oo ) f(x) dx\n";
      cout << "c    is to be approximated by\n";
      cout << "c    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    cout << "c\n";

    for ( i = 0; i < order; i++ )
    {
      cout << "      w(" << i + 1 << ") = " 
           << setprecision(16) << w[i] << "\n";
    }
    cout << "\n";
    for ( i = 0; i < order; i++ )
    {
      cout << "      x(" << i + 1 << ") = " 
           << setprecision(16) << x[i] << "\n";
    }
    cout << "\n";
    for ( i = 0; i < 2; i++ )
    {
      cout << "      r(" << i + 1 << ") = " << r[i] << "\n";
    }
  }
  else if ( output == "F90" )
  {  
    cout << "!\n";
    cout << "!  Weights W, abscissas X and range R\n";
    cout << "!  for a generalized Gauss-Hermite quadrature rule\n";
    cout << "!  ORDER = " << order << "\n";
    cout << "!  ALPHA = " << alpha << "\n";
    cout << "!\n";
    if ( option == 0 )
    {
      cout << "!  OPTION = 0, Standard rule:\n";
	  cout << "!    Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx\n";
	  cout << "!    is to be approximated by\n";
	  cout << "!    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    else
    {
      cout << "!  OPTION = 1, modified rule:\n";
      cout << "!    Integral ( -oo < x < +oo ) f(x) dx\n";
      cout << "!    is to be approximated by\n";
      cout << "!    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    cout << "!\n";

    for ( i = 0; i < order; i++ )
    {
      cout << "  w(" << i + 1 << ") = " 
           << setprecision(16) << w[i] << "\n";
    }
    cout << "\n";
    for ( i = 0; i < order; i++ )
    {
      cout << "  x(" << i + 1 << ") = " 
           << setprecision(16) << x[i] << "\n";
    }
    cout << "\n";
    for ( i = 0; i < 2; i++ )
    {
      cout << "  r(" << i + 1 << ") = " << r[i] << "\n";
    }
  }
  else if ( output == "MAT" )
  {
    cout << "%\n";
    cout << "%  Weights W, abscissas X and range R\n";
    cout << "%  for a generalized Gauss-Hermite quadrature rule\n";
    cout << "%  ORDER = " << order << "\n";
    cout << "%  ALPHA = " << alpha << "\n";
    cout << "%\n";
    if ( option == 0 )
    {
      cout << "%  OPTION = 0, Standard rule:\n";
	  cout << "%   Integral ( -oo < x < +oo ) |x|^ALPHA exp(-x^2) f(x) dx\n";
	  cout << "%   is to be approximated by\n";
	  cout << "%   sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    else
    {
      cout << "%  OPTION = 1, modified rule:\n";
      cout << "%    Integral ( -oo < x < +oo ) f(x) dx\n";
      cout << "%    is to be approximated by\n";
      cout << "%    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).\n";
    }
    cout << "%\n";

    for ( i = 0; i < order; i++ )
    {
      cout << "  w(" << i + 1 << ") = " 
           << setprecision(16) << w[i] << ";\n";
    }
    cout << "\n";
    for ( i = 0; i < order; i++ )
    {
      cout << "  x(" << i + 1 << ") = " 
           << setprecision(16) << x[i] << ";\n";
    }
    cout << "\n";
    for ( i = 0; i < 2; i++ )
    {
      cout << "  r(" << i + 1 << ") = " << r[i] << ";\n";
    }
  } 
  else
  {
    if ( option == 0 )
    {
      output_w = output + "%s_w.txt";
      output_x = output + "%s_x.txt";
      output_r = output + "%s_r.txt";
    }
    else
    {
      output_w = output + "%s_modified_w.txt";
      output_x = output + "%s_modified_x.txt";
      output_r = output + "%s_modified_r.txt";
    }
    cout << "\n";
    cout << "  Creating quadrature files.\n";
    cout << "\n";
    cout << "  Root file name is     \"" << output   << "\".\n";
    cout << "\n";
    cout << "  Weight file will be   \"" << output_w << "\".\n";
    cout << "  Abscissa file will be \"" << output_x << "\".\n";
    cout << "  Region file will be   \"" << output_r << "\".\n";
            
    r8mat_write ( output_w, 1, order, w );
    r8mat_write ( output_x, 1, order, x );
    r8mat_write ( output_r, 1, 2,     r );
  }
  
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void gen_laguerre_compute ( int order, double alpha, double xtab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre quadrature rule.
//
//  Discussion:
//
//    In the simplest case, ALPHA is 0, and we are approximating the
//    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
//    it is easy to modify the rule to approximate the integral from
//    A to +oo as well.
//
//    If ALPHA is nonzero, then there is no simple way to extend the
//    rule to approximate the integral from A to +oo.  The simplest
//    procedures would be to approximate the integral from 0 to A.
//
//    The integration interval is [ A, +oo ) or [ 0, +oo ).
//
//    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x**alpha.
//
//
//    If the integral to approximate is:
//
//        Integral ( A <= X < +oo ) EXP ( - X ) * F(X) dX
//      or
//        Integral ( 0 <= X < +oo ) EXP ( - X ) * X**ALPHA * F(X) dX
//
//    then the quadrature rule is:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
//    or
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
//
//
//    If the integral to approximate is:
//
//        Integral ( A <= X < +oo ) F(X) dX
//      or
//        Integral ( 0 <= X < +oo ) X**ALPHA * F(X) dX
//
//    then the quadrature rule is:
//
//      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) 
//        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
//    or
//      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the quadrature rule to be computed.
//    ORDER must be at least 1.
//
//    Input, double ALPHA, the exponent of the X factor.
//    Set ALPHA = 0.0 for the simplest rule.
//    ALPHA must be nonnegative.
//
//    Output, double XTAB[ORDER], the Gauss-Laguerre abscissas.
//
//    Output, double WEIGHT[ORDER], the Gauss-Laguerre weights.
//
{
  double *b;
  double *c;
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double x;

  b = new double[order];
  c = new double[order];
//
//  Set the recursion coefficients.
//
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( alpha + ( double ) ( 2 * i + 1 ) );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i ) * ( alpha + ( double ) ( i ) );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = r8_gamma ( alpha + 1.0 ) * prod;

  for ( i = 0; i < order; i++ )
  {
//
//  Compute an estimate for the root.
//
    if ( i == 0 )
    {
      x = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) / 
        ( 1.0 + 2.4 * ( double ) ( order ) + 1.8 * alpha );
    }
    else if ( i == 1 )
    {
      x = x + ( 15.0 + 6.25 * alpha ) / 
        ( 1.0 + 0.9 * alpha + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) ) 
        / ( 1.9 * ( double ) ( i - 1 ) );

      r2 = 1.26 * ( double ) ( i - 1 ) * alpha / 
        ( 1.0 + 3.5 * ( double ) ( i - 1 ) );

      ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );

      x = x + ratio * ( x - xtab[i-2] );
    }
//
//  Use iteration to find the root.
//
    gen_laguerre_root ( &x, order, alpha, &dp2, &p1, b, c );
//
//  Set the abscissa and weight.
//
    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;
  }

  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x, 
  int order, double alpha, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_RECUR evaluates a generalized Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Output, double *P2, the value of L(ORDER)(X).
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double X, the point at which polynomials are evaluated.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponent of the X factor in the
//    integrand.
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - alpha - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
//****************************************************************************80

void gen_laguerre_root ( double *x, int order, double alpha, double *dp2, 
  double *p1, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    GEN_LAGUERRE_ROOT improves a root of a generalized Laguerre polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input/output, double *X, the approximate root, which
//    should be improved on output.
//
//    Input, int ORDER, the order of the polynomial to be computed.
//
//    Input, double ALPHA, the exponent of the X factor.
//
//    Output, double *DP2, the value of L'(ORDER)(X).
//
//    Output, double *P1, the value of L(ORDER-1)(X).
//
//    Input, double B[ORDER], C[ORDER], the recursion coefficients.
//
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gen_laguerre_recur ( &p2, dp2, p1, *x, order, alpha, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the 
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  double value;

  value = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + value )  )
  {
    value = value / 2.0;
  }

  value = 2.0 * value;

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for a real argument.
//
//  Discussion:
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
//
//  Coefficients for minimax approximation over (12, INF).
//
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = false;
  fact = one;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}


/*


double f(double alpha, double x)	{
	return 1.0 / (1 + alpha * exp(-x));
}

int main(int argc, char* argv[])	{

	int maxorder = atoi(argv[1]);
	double alpha = atof(argv[2]);

	for (int order = 2; order<= maxorder; order+=2)	{

		double* x = new double[order];
		double* w = new double[order];
		
		gen_hermite_compute(order,0,x,w);

		double total = 0;
		for (int i=0; i<order/2; i++)	{
			total += w[i] * f(alpha,x[i]);
			total += w[order - 1 - i] * f(alpha,x[order - 1 - i]);
		}
		cout.precision(20);
		cout << order << '\t' << total << '\n';
		// cout << order << '\t' << alpha << '\t' << total << '\t' << sqrt(3.141592653589793238)   << '\n';

		delete[] x;
		delete[] w;
	}
	return 0;
}
*/

#endif


