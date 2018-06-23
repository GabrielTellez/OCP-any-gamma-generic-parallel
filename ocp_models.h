#ifndef __OCP_MODELS_H_INCLUDED__
#define __OCP_MODELS_H_INCLUDED__ 

#include <string>
#include <iostream>
#include "ocp_functions.h"
#include <boost/preprocessor/iteration/local.hpp>
#define TOPINDEX_OF_MODELS 6

// include ocp_models_base.h after defining TOPINDEX_OF_MODELS
#include "ocp_models_base.h"

template <class T>
std::string ocp_model<T>::model_names[] = {
  "softdisk",
  "harddisk",
  "sphere_N+1",
  "sphere_N",
  "soft_cylinder_W",
  "manning_Q_Deltainfty",
  "manning_Q_Delta"
};
template <class T>
std::string ocp_model<T>::comment_lines[] = {
  " sqrt{pi rho_b} r      rho(r)/rho_b",
  " sqrt{pi rho_b} r      rho(r)/rho_b",
  " Not implemented",
  " tan(theta/2)          rho(r)/rho_b",
  " y/W                   rho(r)/rho_b",
  " r/R                   rho(r)/rho_b",
  " r/R                   rho(r)/rho_b"
};

// OCP softdisk implementation model_type=0
template <class T>
T ocp_model<T>::f0(unsigned int mu)
{
  return factorial<T>(mu);
}
template <class T>
T ocp_model<T>::f_dens0(unsigned int mu, const T& x)
//  x = sqrt{pi * rho_b} r
{
  // returns (Gamma x^2/2)^{mu}/mu!
  //  std::cout << "Computing f_dens0(mu, x), with mu="<<mu<<", x="<<x<<"\n"; // DEBUG
  return f_dens_softdisk<T>(mu,sqrt(this->gammaover2)*x);  // defined in ocp_functions.h;
}
template <class T>
T ocp_model<T>::pre_f0(const T& x)
{
  return exp(-this->gammaover2*x*x)*this->gammaover2;
}
template <class T>
T ocp_model<T>::free_energy0(const T& Z)
{
  T betaF=-log(Z);
  T Gamma=T(2.0)*T(this->gammaover2);
  T np=T(this->n);
  betaF=betaF-T(3.0)*Gamma*np*np/T(8.0)+(Gamma*np*np/T(4.0))*log(np);
  betaF=betaF-Gamma*np*log(my_pi<T>())/T(4.0);
  betaF=betaF+( (Gamma*np*np/T(4.0)) - np*((Gamma/T(4.0))-T(1.0)) )*log(Gamma/T(2.0));
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model0()
{
  // param1 not used
}

//  OCP Hard disk implementation : model_type=1

template <class T>
T ocp_model<T>::f1(unsigned int mu)
{
  return gammainc<T>(mu+1,T(this->gammaover2*this->n));
}
template <class T>
T ocp_model<T>::f_dens1(unsigned int mu, const T& x)
//  x = sqrt{pi * rho_b} r
{
  // returns (Gamma x^2/2)^{mu}/gammainc(mu+1, Gamma * n/2)
  if ( x*x > T(this->n) ) {
    return 0.0;
  } else {
    return my_pow<T>(T(this->gammaover2)*x*x,mu)/gammainc<T>(mu+1, T(this->gammaover2*this->n));
  }
}
template <class T>
T ocp_model<T>::pre_f1(const T& x)
{
  return exp(-this->gammaover2*x*x)*this->gammaover2;
}
template <class T>
T ocp_model<T>::free_energy1(const T& Z)
{
  T betaF=-log(Z);
  T Gamma=T(2.0)*T(this->gammaover2);
  T np=T(this->n);
  betaF=betaF-T(3.0)*Gamma*np*np/T(8.0)+(Gamma*np*np/T(4.0))*log(np);
  betaF=betaF-Gamma*np*log(my_pi<T>())/T(4.0);
  betaF=betaF+( (Gamma*np*np/T(4.0)) - np*((Gamma/T(4.0))-T(1.0)) )*log(Gamma/T(2.0));
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model1()
{
  // param1 not used
}

// OCP sphere N+1 particles implementation model_type=2
// Full partition function for N+1 particles is
// Q= Z * exp(Gamma(N+1)^2/4) pi^(N+1) * 1/((N+1)*(2R)^(2(N+1)(Gamma/4-1)))
//      
// Program computes Z
//
template <class T>
T ocp_model<T>::f2(unsigned int mu)
{
  //  std::cout << mu << " " << this->gammaover2 << " " << this->n <<"\n" ; //DEBUG
  T tmp=factorial<T>(mu+(this->gammaover2))*factorial<T>(((this->n)-1)*(this->gammaover2)-mu);
  tmp/=factorial<T>((this->n)*(this->gammaover2)+1);
  return tmp;
}
template <class T>
T ocp_model<T>::f_dens2(unsigned int mu, const T& x)
//  not implemented yet
{
  return 0.0;
}
template <class T>
T ocp_model<T>::pre_f2(const T& x)
//  not implemented yet
{
  return 0.0;
}
template <class T>
T ocp_model<T>::free_energy2(const T& Z)
{
  T betaF=-log(Z);
  T np=T(this->n + 1);
  betaF=betaF - T(this->gammaover2)*np*np/T(2.0);
  betaF=betaF - np*T(this->gammaover2)*log(my_pi<T>())/T(2.0);
  betaF=betaF + np*((T(this->gammaover2)/T(2.0))-T(1.0))*log(np);
  betaF=betaF+log(np);
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model2()
{
  // param1 not used
}


// OCP sphere N particles implementation model_type=3
// Full partition function for N particles is
// Q= Z * exp(Gamma N^2/4) pi^N /((2R)^(2N(Gamma/4-1)))
//      
// Program computes Z
//
template <class T>
T ocp_model<T>::f3(unsigned int mu)
{
  //  std::cout << mu << " " << this->gammaover2 << " " << this->n <<"\n" ; //DEBUG
  T tmp=factorial<T>(mu)*factorial<T>(((this->n)-1)*(this->gammaover2)-mu);
  tmp/=factorial<T>( ((this->n)-1)*(this->gammaover2)+1);
  return tmp;
}
template <class T>
T ocp_model<T>::f_dens3(unsigned int mu, const T& x)
{
  return my_pow<T>(x*x,mu)/f3(mu);
}
template <class T>
T ocp_model<T>::pre_f3(const T& x)
{
  T tmp = 1.0;
  tmp/=T(this->n)*(my_pow<T>(1+x*x,(this->n-1)*this->gammaover2));
  return tmp;
}
template <class T>
T ocp_model<T>::free_energy3(const T& Z)
// Full partition function for N particles is
// Q= Z * exp(Gamma N^2/4) pi^N /((2R)^(2N(Gamma/4-1)))
// beta F = - log(Q)
{
  T betaF=-log(Z);
  T np=T(this->n);
  betaF=betaF - T(this->gammaover2)*np*np/T(2.0);
  betaF=betaF - np*T(this->gammaover2)*log(my_pi<T>())/T(2.0);
  betaF=betaF + np*((T(this->gammaover2)/T(2.0))-T(1.0))*log(np);
  //  betaF=betaF+log(np);
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model3()
{
  // param1 not used
}


// OCP soft cylinder
// Partition function (reduced, see Forrester draft) for N particles is
// Q = Z * N! (W/Gamma)^{N/2} / rho_b^N
//      
// Program computes Z
//
template <class T>
T ocp_model<T>::f4(unsigned int mu)
{
  T tmp=exp(my_pi<T>() *2.0*T(this->gammaover2 )*(((T(this->n)/2.0)-((T(mu)/T(this->gammaover2))+0.5))*((T(this->n)/2.0)-((T(mu)/T(this->gammaover2))+0.5))) /this->param1);
  return tmp;
}
template <class T>
T ocp_model<T>::f_dens4(unsigned int mu, const T& x)
// x = y/W = y_tilde/param1 
{
  return exp(-4.0 * my_pi<T>()*x*T(mu) ) /f4(mu);
}
template <class T>
T ocp_model<T>::pre_f4(const T& x)
{
  T tmp = exp(- my_pi<T>() * 2.0 * T(this->gammaover2) * x * x / T(this->param1) );
  tmp = tmp * exp(2.0*T(this->gammaover2)*my_pi<T>()*T(this->n-1)* x);

  return tmp;
}
template <class T>
T ocp_model<T>::free_energy4(const T& Z)
{
  T betaF = -log(Z);
  betaF = betaF-T(this->n)*T(this->gammaover2)*log(T(2.0)*my_pi<T>());
  betaF = betaF+T(this->n)*T(this->gammaover2)*log(this->param1)/T(2.0);
  betaF = betaF+T(this->n*this->n*this->n)*T(this->gammaover2)*my_pi<T>()/(T(6.0)*this->param1);
  betaF = betaF - T(this->n)*log(this->param1)/T(2.0);
  betaF = betaF + T(this->n)*log(T(this->gammaover2)*T(2.0))/T(2.0);
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model4()
{
  // param1 = W_tilde = cylinder perimeter in reduced units.
  if (this->param1 <= 0.0 ) {
    std::cerr << "param1 = W_tilde = perimeter in reduced units = " << this->param1 << " cannot be negative.\n";
    exit(-1);
  }
}
// OCP Manning model: Charged disk with charge Q and N counterions charge -1
// in a cell with D->infty
// Q = xi/Xi = 2 xi /Gamma = param1
// Program computes Z
// Z = sum_mu (c_mu^2/prod mi!) prod (1/(Q*Gamma/2-mu_l-1))
//
template <class T>
T ocp_model<T>::f5(unsigned int mu)
{
  T tmp=T(this->gammaover2)*T(this->param1)-T(mu+1);
  tmp = 1.0/tmp;
  return tmp;
}
template <class T>
T ocp_model<T>::f_dens5(unsigned int mu, const T& x)
// x = r/R 
{
  return my_pow<T>(x*x,mu)/f5(mu);
}
template <class T>
T ocp_model<T>::pre_f5(const T& x)
{
  T tmp = 1.0/my_pi<T>();
  tmp = tmp*exp((-T(2.0)*this->param1*T(this->gammaover2))*log(x));
  return tmp;
}
template <class T>
T ocp_model<T>::free_energy5(const T& Z)
{
  // R=1, L=1  
  T betaF = -log(Z)-T(this->n)*log(my_pi<T>());
  //  betaF=betaF+(this->gammaover2)((this->n)-(this->param1))*(this->param2);
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model5()
{
  // param1 = Q = colloid charge
  if (this->param1 - T(this->n) + T(1.0) - T(1.0/this->gammaover2) < 0.0 ) {
    std::cerr << "Warning:\n";
    std::cerr << "param1 = Q = colloid charge = " << this->param1 << " .\n";
    std::cerr << "Q must satisfy Q+1-N-(2/Gamma) = ";
    std::cerr << this->param1 - T(this->n) + T(1.0) - T(1.0/this->gammaover2) << " > 0\nProceeding anyway...";
  }
  if (this->param1 - T(this->n) + T(1.0) - T(1.0/this->gammaover2) < 0.0 ) {
    std::cerr << "param1 = Q = colloid charge = " << this->param1 << " .\n";
    std::cerr << "Q must satisfy Q+1-N-(2/Gamma) = ";
    std::cerr << this->param1 - T(this->n) + T(1.0) - T(1.0/this->gammaover2) << " != 0\n";
    exit(-1);
  }
}
// OCP Manning model: Charged disk with charge Q and N counterions charge -1
// in a cell of radius D finite
// param1 = Q =xi /Xi
// param2 = Delta = log(D/R)  ; with  R=1
// Program computes Z
// Z = sum_mu (c_mu^2/prod mi!) prod ( (1-(D/R)^2(mu_l+1-Gamma *Q/2))/(Q*Gamma/2-mu_l-1))
//
template <class T>
T ocp_model<T>::f6(unsigned int mu)
{
  if ( T(mu+1)-T(this->gammaover2*this->param1) == 0 ) {
    return 2*(this->param2);
  } else {
    T tmp=T(this->gammaover2)*T(this->param1)-T(mu+1);
    tmp = (1.0-exp(2.0*(T(mu+1)-T(this->gammaover2)*this->param1)*(this->param2)))/tmp;
    return tmp;
  }
}
template <class T>
T ocp_model<T>::f_dens6(unsigned int mu, const T& x)
// x = r/R 
{
  return my_pow<T>(x*x,mu)/f6(mu);
}
template <class T>
T ocp_model<T>::pre_f6(const T& x)
{
  T tmp = 1.0/my_pi<T>();
  tmp = tmp*exp((-T(2.0)*this->param1*T(this->gammaover2))*log(x));
  return tmp;
}
template <class T>
T ocp_model<T>::free_energy6(const T& Z)
{
  // R=1, L=1  
  T betaF = -log(Z)-T(this->n)*log(my_pi<T>());
  betaF=betaF+T(this->gammaover2)*(T(this->n)-(this->param1))*(T(this->n)-(this->param1))*(this->param2);
  return betaF;
}
template <class T>
void ocp_model<T>::validate_model6()
{
  // param1 = Q = colloid charge
  if (this->param2 <= 0.0) {
    std::cerr << "p2 = log(D/R) = " << this->param2 << " = Wigner cell size must be greater than 0.\n";
    exit(-1);
  }
}




#endif




