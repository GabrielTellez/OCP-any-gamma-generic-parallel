#ifndef __HARDDISK_FUNCTIONS_H_INCLUDED__
#define __HARDDISK_FUNCTIONS_H_INCLUDED__ 

#include "softdisk_functions.h"

template <class T>
T f_harddisk(unsigned int mu, long int n_gam)
// gamma(mu+1,n_gam)
{
  if (mu==0) {
    return T(1)-exp(-T(n_gam)); 
  } else {
    return T(mu)*f_harddisk<T>(mu-1,n_gam)-exp(-T(n_gam))*my_pow(T(n_gam),mu);
  }  
}

// To do....
template <class T>
T f_dens_harddisk(unsigned int mu, const T& x)  
//  x = sqrt{gammaover2* pi * rho_b} r
{
  return my_pow(x,2*mu);
}

template <class T>
T pre_f_harddisk(const T& x)
{
  return exp(-x*x);
}

#endif
