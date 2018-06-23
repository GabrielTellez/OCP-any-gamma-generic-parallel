#ifndef __OCP_MODELS_BASE_H_INCLUDED__
#define __OCP_MODELS_BASE_H_INCLUDED__ 

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <boost/preprocessor/iteration/local.hpp>

#include "process_options.h"

template <class T>
class ocp_model_basic
{
 public:
  //  ocp_model_basic(long int n_i, int g_i);
  ocp_model_basic(process_options& op);
  std::string get_name() { return name; }
  std::string get_comment_line() { return comment_line; }
 protected:
  long int n;  // number of particles
  int gammaover2;   // coupling Gamma/2
  T param1;  // parameter for some models (softdisk W, Manning 2D Q)
  T param2;  // parameter for some models (Manning 2D D cell radius)
  unsigned int model_type;
  std::string name; // name of the model
  std::string comment_line;   // for density files
};

template <class T>
ocp_model_basic<T>::ocp_model_basic(process_options& op)
: n(op.n), gammaover2(op.gammaover2), param1(op.param1),
  param2(op.param2), model_type(op.model)
{}

template <class T>
class ocp_model : public ocp_model_basic<T>
{
 public:
  ocp_model(process_options& op);
  std::string get_name(); 
  T f(unsigned int);  // function that goes in the product over partition elements
  T f_dens(unsigned int, const T&);  // function to compute the density that goes in the sum over partition elements
  T pre_f(const T&); // prefactor function to normalize densities
  T free_energy(const T&); // computes the free energy from the partition function
  void validate_model(); // validates if param1 is OK for model type
  static void print_models(); // prints all available models
 protected:
  static std::string model_names[TOPINDEX_OF_MODELS+1];  
  static std::string comment_lines[TOPINDEX_OF_MODELS+1]; 
  typedef T (ocp_model::*type_f)(unsigned int);
  type_f jumptable_f[TOPINDEX_OF_MODELS+1];
  typedef T (ocp_model::*type_f_dens)(unsigned int, const T&);
  type_f_dens jumptable_f_dens[TOPINDEX_OF_MODELS+1];
  typedef T (ocp_model::*type_pre_f)(const T&);
  type_pre_f jumptable_pre_f[TOPINDEX_OF_MODELS+1];
  typedef T (ocp_model::*type_free_energy)(const T&);
  type_pre_f jumptable_free_energy[TOPINDEX_OF_MODELS+1];
  typedef void (ocp_model::*type_validate_model)();
  type_validate_model jumptable_validate_model[TOPINDEX_OF_MODELS+1];
  #define BOOST_PP_LOCAL_MACRO(n) \
  T f##n(unsigned int);		  \
  T f_dens## n(unsigned int, const T&);  \
  T pre_f##n(const T&);                 \
  T free_energy##n(const T&);                 \
  void validate_model##n();                 \
   /**/
  #define BOOST_PP_LOCAL_LIMITS (0, TOPINDEX_OF_MODELS)
  #include BOOST_PP_LOCAL_ITERATE()

};

template <class T>
ocp_model<T>::ocp_model(process_options& op)
: ocp_model_basic<T>(op)
{
  this->name=model_names[this->model_type];
  this->comment_line=comment_lines[this->model_type];
#define BOOST_PP_LOCAL_MACRO(n)	          \
  jumptable_f[n]=&ocp_model<T>::f##n ;	  \
  jumptable_f_dens[n]=&ocp_model<T>::f_dens##n ;	  \
  jumptable_pre_f[n]=&ocp_model<T>::pre_f##n ;		  \
  jumptable_free_energy[n]=&ocp_model<T>::free_energy##n ;		  \
  jumptable_validate_model[n]=&ocp_model<T>::validate_model##n ;  \
    /**/
#define BOOST_PP_LOCAL_LIMITS (0, TOPINDEX_OF_MODELS)
#include BOOST_PP_LOCAL_ITERATE()

  validate_model();
} // Constructor ocp_model
 
template <class T>
std::string ocp_model<T>::get_name()
{
  std::stringstream myname;
  //  myname <<  setiosflags(std::ios::scientific) << std::setprecision(12) ;
  myname << this->name; 
  if (this->param1 !=0) {
    myname << this->param1 ;
    if(this->param2 != 0) {
    myname << "_" << this->param2 ;
    }
  } else {
    if(this->param2 != 0) {
    myname << this->param1 ;
    myname << "_" << this->param2 ;
    }
  }
  return myname.str();
}

template <class T>
T ocp_model<T>::f(unsigned int mu)
{
  return (this->*(jumptable_f[this->model_type]))(mu);
}
template <class T>
T ocp_model<T>::f_dens(unsigned int mu, const T& x)
{
  return (this->*(jumptable_f_dens[this->model_type]))(mu,x);
}
template <class T>
T ocp_model<T>::pre_f(const T& x)
{
  return (this->*(jumptable_pre_f[this->model_type]))(x);
}
template <class T>
T ocp_model<T>::free_energy(const T& Z)
{
  return (this->*(jumptable_free_energy[this->model_type]))(Z);
}
template <class T>
void ocp_model<T>::validate_model()
{
  (this->*(jumptable_validate_model[this->model_type]))();
}

template <class T>
void ocp_model<T>::print_models()
{
  std::cout << "Available models are:" << std::endl;
  for( unsigned int i=0; i<= TOPINDEX_OF_MODELS ; i++) {
    std::cout << i << " : " << model_names[i] << std::endl;
  }
}



#endif
