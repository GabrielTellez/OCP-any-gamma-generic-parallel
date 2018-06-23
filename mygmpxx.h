#ifndef __MYGMPXX_H_INCLUDED__
#define __MYCMPXX_H_INCLUDED__ 

#include <sstream>
#include <gmpxx.h>

// define a convertion from mpf_class to long double

/* This very slow 
mpf_class ld_to_mpf_class( long double x )
{
  std::stringstream str;
  str << x;
  return mpf_class(str.str());
}
*/

template <class T>
T ld_to_mp(long double x)
{
  return static_cast<T>(x);
}

// Let's hope (long double) coefz holds in a double
template <>
mpf_class ld_to_mp(long double x)
{
  return static_cast<mpf_class>(double(x));
}


#endif
