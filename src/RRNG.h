// YOU MUST ALWAYS CALL GetRNGSeed() and PutRNGSeed() WHEN USING THESE FUNCTIONS!!!

#ifndef __BASICRNG__
#define __BASICRNG__

#include "R.h"
#include "Rmath.h"
#include "Matrix.h"

class BasicRNG {

 public:

  // Random variates.
  inline double unif  ();                             // Uniform
  inline double expon_mean(double mean);                  // Exponential
  inline double expon_rate(double rate);                  // Exponential
  inline double chisq (double df);                    // Chisq
  inline double norm  (double sd);                    // Normal
  inline double norm  (double mean , double sd);      // Normal
  inline double gamma_scale (double shape, double scale); // Gamma_Scale
  inline double gamma_rate  (double shape, double rate);  // Gamma_Rate
  inline double igamma(double shape, double scale);   // Inv-Gamma
  inline double flat  (double a=0  , double b=1  );   // Flat

  inline int bern  (double p);                     // Bernoulli

  // CDF
  static inline double p_norm (double x, int use_log=0);

}; // BasicRNG

//////////////////////////////////////////////////////////////////////
		      // R Random Variates //
//////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------
// Distributions with one parameter.

#define ONEP(NAME, CALL, P1)			\
  inline double BasicRNG::NAME(double P1)	\
  {						\
    return CALL (P1);				\
  }						\

ONEP(expon_mean, rexp  , mean)
ONEP(chisq     , rchisq, df  )

#undef ONEP

//--------------------------------------------------------------------
// Distributions with two parameters.

#define TWOP(NAME, CALL, P1, P2)			\
  inline double BasicRNG::NAME(double P1, double P2)	\
  {							\
    return CALL (P1, P2);				\
  }							\

TWOP(gamma_scale, rgamma, shape, scale)
TWOP(norm      , rnorm , mean , sd  )
TWOP(flat      , runif , a    , b   )

// x ~ Gamma(shape=a, scale=b)
// x ~ x^{a-1} exp(x / b).

#undef TWOP

//--------------------------------------------------------------------
			    // Uniform //

inline double BasicRNG::unif()
{
  return unif_rand();
} // unif

//--------------------------------------------------------------------
			  // Exponential //
inline double BasicRNG::expon_rate(double rate)
{
  return expon_mean(1.0 / rate);
}

//--------------------------------------------------------------------
			    // Normal //

inline double BasicRNG::norm(double sd)
{
  return rnorm(0, sd);
} // norm

//--------------------------------------------------------------------
			   // gamma_rate //

inline double BasicRNG::gamma_rate(double shape, double rate)
{
  return gamma_scale(shape, 1.0 / rate);
}

//--------------------------------------------------------------------
			   // Inv-Gamma //

// a = shape, b = scale
// x ~ IG(shape, scale) ~ x^{-a-1} exp(b / x).
// => 1/x ~ Ga(shape, scale*=1/scale).

inline double BasicRNG::igamma(double shape, double scale)
{
  return 1.0/rgamma(shape, 1.0 / scale);
} // igamma

////////////////////////////////////////////////////////////////////////////////

inline double BasicRNG::p_norm(double x, int use_log)
{
  return pnorm(x, 0.0, 1.0, 1, use_log);
}

////////////////////////////////////////////////////////////////////////////////

#endif
