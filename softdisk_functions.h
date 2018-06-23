#ifndef __SOFTDISK_FUNCTIONS_H_INCLUDED__
#define __SOFTDISK_FUNCTIONS_H_INCLUDED__ 

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
T f_softdisk(unsigned int mu)
// = factorial(mu) = Gamma(mu+1)
{
  //  return gamma(x+1.0);  // Broken on arprec-2.2.14  (painfully slow); OK for 2.2.6
  
  if (mu==0) {
    return T(1); 
  } else {
    return T(mu)*f_softdisk<T>(mu-1);
  }  
}

template <class T>
T f_dens_softdisk(unsigned int mu, const T& x)  
//  x = sqrt{gammaover2* pi * rho_b} r
{
  // returns x^{2 mu}/mu!
  // recursive: x^{2 mu}/mu! = (x^2/mu)* (x^2)^{mu-1}/(mu-1)!
  if (mu==0) {
    return 1;
  } else {
    return x*x*f_dens_softdisk(mu-1,x)/T(mu);
  }
}

template <class T>
T f_dens_softdisk_g(unsigned int mu, const T& x, int gammaover2)  
//  x = sqrt{pi * rho_b} r
{
  // returns (Gamma x/2)^{2 mu}/mu!
  return f_dens_softdisk(mu,sqrt(gammaover2)*x);
}

template <class T>
T pre_f_softdisk(const T& x)
{
  return exp(-x*x);
}

template <class T>
T pre_f_softdisk_g(const T& x, int gammaover2)
{
  return exp(-gammaover2*x*x);
}


#endif
