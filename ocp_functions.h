#ifndef __OCP_FUNCTIONS_H_INCLUDED__
#define __OCP_FUNCTIONS_H_INCLUDED__ 
#include <arprec/mp_real.h>
#include <math.h>

template <class T>
T my_pi();


template <class T>
T factorial(unsigned int mu)
{
  //  std::cout << "factorial arg = " << mu << "\n" ; //DEBUG
  if (mu==0) {
    return T(1); 
  } else {
    return T(mu)*factorial<T>(mu-1);
  }  
}

template <class T>
T f_dens_softdisk(unsigned int mu, const T& x)  
//  x = sqrt{gammaover2* pi * rho_b} r
{
  // returns x^{2 mu}/mu!
  // recursive: x^{2 mu}/mu! = (x^2/mu)* (x^2)^{mu-1}/(mu-1)!
  //  std::cout << "computing f_dens_softdisk(mu,x) with mu= "<< mu << ", x=" << x << "\n";
  if (mu==0) {
    return T(1.0);
  } else {
    if (x==0) {
      return T(0.0);
    } else {
      return x*x*f_dens_softdisk<T>(mu-1,x)/T(mu);
    }
  }
}

template <class T>
T my_pow(const T& x, unsigned int mu)
{
  if(mu==0) {
    return T(1);
  } else {
    return x*my_pow(x, mu-1);
  }

}

template <class T>
T gammainc(unsigned int mu, const T& z)
// return incomplete gamma function gamma(mu,z)
{
  if (mu==0) {
    std::cerr << "gammainc: error: first argument mu cannot be 0\n";
    exit(1);
  }
  T expz=exp(-z);
  if ( mu==1 ) {
    return 1.0-expz;
  } else {
    return T(mu-1)*gammainc(mu-1,z)-expz*my_pow(z,mu-1);
  }
}

template <class T>
T HarmonicNumber(unsigned int mu)
// returns the harmonic number mu, used in the computation of the sphere internal energy
{
  if (mu==0) {
    std::cerr << "Harmonic number: error: first argument mu cannot be 0\n";
    exit(1);
  }
  if ( mu==1 ) {
    return T(1.0);
  } else {
    return 1.0/T(mu)+HarmonicNumber<T>(mu-1);
  }
}

#endif
