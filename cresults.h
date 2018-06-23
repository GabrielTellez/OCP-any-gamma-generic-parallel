#ifndef __CRESULTS_H_INCLUDED__
#define __CRESULTS_H_INCLUDED__ 

#include <string>
#include <sstream>
#include <vector>
#include <arprec/mp_real.h>
#include "partition.h"
#include "process_options.h"

// #include <boost/ptr_container/ptr_vector.hpp>


template <class T>
class Domain
{
protected:
  T xmin;
  T xmax;
  size_t xsize;
public:
  Domain<T> (); // Default constructor with xmin=xmax=0, xsize=2
  Domain<T> (const T& xmin_i, const T& xmax_i, const size_t  xsize_i);
  Domain<T> (const std::string s);  // constructor from string s="min,max,size"
  Domain<T> (const Domain<T> & dom);
  T min() const {return xmin; }
  T max() const { return xmax;}
  size_t size() const { return xsize; }
  T step() const { return (xmax-xmin)/(int(xsize-1)); }
  // Gives the point x corresponding to the index i
  T x(size_t i) const; 
  // compare if two domains are equal
  bool operator==(const Domain<T>& dom) const;
  void set_xmin(T x) { xmin=x; }
  void set_xmax(T x) { xmax=x; }
  void set_size( size_t s) { xsize=s;   
    if(size() <2 ) {
      std::cerr << "Domain size has to be >=2, setting it to 2\n";
      xsize=2;
    } 
  }
  void set(std::string); // set the domain with a string formated as: "xmin,xmax,xsize"
  void print_domain(std::ostream& os, std::string sep="_");
 // Gives the index i corresponding to x
  //  size_t index(const T& x) const { return size_t((x-xmin)/(this->step()));}
};
template <class T>
std::istream& operator>> (std::istream& in, Domain<T>& dom );
template <class T>
std::ostream& operator<< (std::ostream& os, Domain<T>& dom );



template <class T>
class DomainId : public Domain<T>
{
protected:
  std::string name;
public:
  DomainId<T>(const Domain<T>& dom_i, const std::string name_i="");
  DomainId<T>(const T& xmin_i, const T& xmax_i, const size_t xsize_i, const std::string name_i="");
  DomainId<T>(const std::string s, std::string name_i="");  // constructor from string s="min,max,size"
  DomainId<T> (const DomainId<T> & dom_i);
  std::string get_name() const { return name ; }
  void set_name(std::string name_i);
};



// Precomputed function of two arguments an integer and a T  
template <class T>
class function 
{
 protected:
  long int n; // maximum integer for which the function will be evaluated
 public:
  DomainId<T> dom;
  size_t index_x;  // value of index_x to be bound
 private:
  std::vector<std::vector<T> > table; // precomputed table of values table[mu][i] with i index corresponding to x
 public:
  template <class F> function<T>(long int n_i , DomainId<T> dom_i, F fnc);
  //  T operator() (long int mu, T x) const; // return function value with x given
  T operator() (long int mu, size_t index_x_i) const; // return function value with index given
  T operator() (long int mu) const;
};

template <class T>
template <class F> 
function<T>::function(long int n_i , DomainId<T> dom_i, F fnc)
: n(n_i), dom(dom_i), table(n+1, std::vector<T>(dom.size(),T(0.0)))
{
  // initialize the table of precomputed values
  for(unsigned int mu=0; mu<=size_t(n); mu++) {
    for(unsigned int i=0; i< dom.size(); i++) {
      table[mu][i]=fnc(mu, dom.x(i));
    }
  }
}



// Precomputed function of an integer
template <class T>
class unifunction 
{
 protected:
  std::vector<T> table;
 public:
  template <class F> unifunction<T>(long int n_i, F);
  unifunction<T>(const unifunction<T>&);
  T operator() (long int mu) const;
};

template <class T>
template <class F>
unifunction<T>::unifunction(long int n_i, F fnc)
  : table(n_i+1,0.0)
{
  // initialize the table of precomputed values
  for(unsigned int mu=0; mu<table.size(); mu++) {
    // debug    std::cout  << "Computing f(mu) for mu = " << mu  <<"\n";
    table[mu]=fnc(mu);
  }
}


// holds a set of pointers to functions F=function<T>

template <class T>
class function_set
{
 protected:
  std::vector<function<T>* > table;
 public:
  function_set();
  ~function_set();
  template <class Func> function_set(std::vector<std::string> domains, Func, std::string filenamebase, int gammaover2, long int n, std::string mp_type="");
  function<T>* operator[] (size_t);
  size_t size() const;
};
template <class T>
function_set<T>::function_set() 
{}
template <class T>
template <class Func>
function_set<T>::function_set(std::vector<std::string> domains, Func f, std::string filenamebase, int gammaover2, long int n, std::string mp_type)
{
  for(size_t i = 0 ; i <  domains.size() ; i++ ) {
    std::stringstream filename;
    filename << filenamebase << "_r_" << domains.at(i)
	     << "_G" << 2*gammaover2 << "_n" << n ;
    if (mp_type != "") {
      filename << "_" << mp_type;
    }
    table.push_back(new function<T>((n-1)*gammaover2, 
				DomainId<T>(domains.at(i), filename.str()), f));
    // DEBUG    std::cout << "created function " << i << "\n" << std::flush ;
  }
  //  DEBUG std::cout << "created " << size()<<  "functions.\n" << std::flush ;
}

template <class T>
function_set<T>::~function_set()
{
  //  std::cout << "~function_set()\n" ;  // DEBUG
  for (size_t i=0; i< size() ; i++) {
    delete table.at(i);
  }
}

template <class T>
function<T>* function_set<T>::operator[](size_t i) 
{
  return table.at(i);
}
template <class T>
size_t function_set<T>::size() const
{
  return table.size();
}






// class to hold the computed density data
template <class T>
class Cdensity_result : public DomainId<T>
{
 public:
  std::vector<T> density;
 public:
  Cdensity_result<T>(const DomainId<T>& dom_i, std::string comment_i="");
  Cdensity_result<T>(const Cdensity_result<T>& a);
  T operator[](size_t i);
  bool compare_type(const Cdensity_result<T>& a) const;
  Cdensity_result<T> operator+ (const Cdensity_result<T>& a) const;
  Cdensity_result<T>& operator+= (const Cdensity_result<T>& a);
  Cdensity_result<T> operator/(const T& z) const;
  Cdensity_result<T> normalize(const T& z);
  // Multiplies the density by a function of x
  Cdensity_result<T> operator*=(T (*fnc)(const T&)); 
  template <class F> Cdensity_result<T> normalize(F, const T& z);
  void print(std::string filename, int prec);
  void print(std::ostream& os, int prec);
  // updates the density with data from partition part
  void update(const T& prefactor, cpartition& part, function<T>* fnc_dens);
  void clear();
 protected:
  std::string comment;  // comment line to print in the first line of the archive when saving (print method)
};

template <class T>
template <class F>
Cdensity_result<T> Cdensity_result<T>::normalize(F fnc, const T& z)
{
  for (size_t x_index = 0 ; x_index < density.size() ; x_index++ ) {
    (density.at(x_index))*=fnc(this->x(x_index))/z;
  }
  return *this;
}



// Class that holds a density and a pointer to the function used to compute the density
template <class T>
class Cdensity : public Cdensity_result<T>
{

 protected:
  function<T>* fnc_dens; // function pointer (precomputed) used in the computation of the density (in sum over partitions)
  //  typedef T(*F)(const T&);
  //  F fnc_pref;
  //  T (*fnc_pref)(const T&); // function pointer to normalize the final density (simple prefactor: ie exp(-x*x) for soft disk)
 public:
  Cdensity<T>(const DomainId<T>& dom_i, function<T>* func_dens_i, std::string com="");
  Cdensity<T>(const Cdensity<T>& a); // copy constructor
  Cdensity<T>(const Cdensity<T>* dens_ptr); // copy only the functions pointers, but the density is initialized to 0.0
  void update(const T& prefactor, cpartition& part);
  function<T>* get_fnc_ptr() { return fnc_dens; }
};

template <class T>
  Cdensity<T>::Cdensity(const DomainId<T>& dom_i, function<T>* f_i, std::string com) 
  : Cdensity_result<T> (dom_i, com), fnc_dens(f_i)
{ }

// Class Cdensity (END)


//Class that holds a set of densities
template <class T>
  class Cdensity_set : public std::vector<Cdensity<T> >
{
 public:
  Cdensity_set<T> (); // create empty set
  Cdensity_set(const Cdensity_set<T>& a); // copy constructor: copy all data from a
  Cdensity_set(Cdensity_set<T>* set_ptr); // copy function pointers from set_ptr but densities are cleared to zero
  Cdensity_set<T> operator+ (const Cdensity_set<T>& a) const;
  Cdensity_set<T>& operator+= (const Cdensity_set<T>& a);
  template <class F> Cdensity_set<T> normalize(F, const T& z);
  void print(std::string filename, int prec);
  void update(const T& prefactor, cpartition& part);
  void initialize_new_density(function<T>* f_i, std::string com="");
  void initialize_density_set(function_set<T>& , std::string com="");
};

template <class T>
template <class F>
Cdensity_set<T> Cdensity_set<T>::normalize(F f, const T& z)
{
  for (size_t i=0; i<this->size(); i++) {
    (this->at(i)).normalize(f, z);
  }
  return *this;
}



// Class Model
template <class T>
class Model
{
 public:
  long int n; // number of particles
  int gammaover2; // coupling constant Gamma/2
  long double checksum;
  long double Z_ld; // partition_function in long double precision
  T Z_mp; // partition function in arbitrary precision
  // To add: densities and correlations

  Model(long int n_i, int gam_i, unifunction<long double>* fnc_ld_i, unifunction<T>* fnc_mp_i);
// fnc = unifunction used to compute the partition function
  Model(Model * mod);
  // copy model data from mod but with empty results = 0.0

  Model& operator+= (const Model& a);
  // update results with data from partition part
  Model update(long double coefz, cpartition& part);
  Model operator+ (const Model& a) const;

 protected:
  unifunction<long double>* fnc_ld;
  unifunction<T>* fnc_mp;

};


// Basic Model that computes the checksum and the partition function
template <class T>
class ModelBase
{
 public:
  long int n; // number of particles
  int gammaover2; // coupling constant Gamma/2
  T count; // counts the total number of partitions
  T checksum;
  T Z_mp; // partition function in arbitrary precision
  T free_energy; // free energy
  int max_prec;  // precision used for computations

  ModelBase(long int n_i, int gam_i, unifunction<T>* fnc_mp_i,int max_prec_i=14);
// fnc = unifunction used to compute the partition function

  ModelBase(const ModelBase& mod);
  // copy all model data from mod

  ModelBase(ModelBase * mod);
  // copy model data from mod but with empty results = 0.0

  ModelBase& operator+= (const ModelBase& a);

  // update results with data from partition part, and returns the prefactor
  // used to compute the partition function, 
  // so it can be used in inherited method to compute the density
  T update(cpartition& part);
  ModelBase operator+ (const ModelBase& a) const;

  void print_report(std::ostream& os, process_options& op, std::string mp_type="", std::string model_name="");
 protected:
  unifunction<T>* fnc_mp;

};

// Model that computes the checksum, partition function (from ModelBase) and one density
template <class T>
class Model_unidensity : public ModelBase<T>
{
 public:
  Cdensity_result<T> density;
  Model_unidensity(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, function<T>* fnc_dens_i, int max_prec_i);
  Model_unidensity(Model_unidensity * mod);
  Model_unidensity& operator+= (const Model_unidensity& a);
  // update results with data from partition part
  T update(cpartition& part);
  Model_unidensity operator+ (const Model_unidensity& a) const;
  // divide the density by the partition function and multiply by fnc(x) * factor
  void normalize_density(T (*fnc)(const T&)); 
  void normalize_density(T (*fnc)(const T&), const T& factor = 1);
 protected:
  function<T>* fnc_dens;  // To compute the density at x

};

// Model that computes the checksum, partition function and a set of densities
// Class Density_Set should have defined operator+=, +, update, print, normalize
// for example Density_Set = Cdensity_set<T>
//
template <class T, class Density_Set = Cdensity_set<T> >
class Model_density : public ModelBase<T>
{
 public:
 Density_Set density;
 Model_density(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, int max_prec_i);
 Model_density(const Model_density& mod); // copy all data from mod
 Model_density(Model_density * mod); // copy function pointers from *mod, but densities=0.0
 Model_density<T, Density_Set>& operator+= (const Model_density<T, Density_Set>& a);
 Model_density<T, Density_Set> operator+ (const Model_density<T, Density_Set>& a) const;
 // update results with data from partition part
 T update(cpartition& part);
 // divide the density by the partition function and multiply by fnc(x) * factor
 template <class F> void normalize_density(F, const T& factor=1);
 void initialize_new_density(function<T>* f_i, std::string com=""); 
 void initialize_density_set(function_set<T>& , std::string com="");
 void output_report(process_options& op, std::string mp_type, std::string model_name);
};

template <class T, class Density_Set>
template <class F>
void Model_density<T, Density_Set>::normalize_density(F fnc, const T& factor)
{
  density.normalize(fnc, (this->Z_mp)/factor);
}


// Class that holds the coefficients of the polynomial expansion of the density and the pair correlation function
template <class T>
class Density_coefs
{
 protected:
  int gammaover2;  // Gamma/2
  long int n;  // number of particles N
  unifunction<T>* fnc; // unifunction pointer (precomputed) used in the computation of the density (same as the one for the partition function
 public:
  std::vector<T> density;   // coefficients of the density
  std::vector<T> correlation;  // coefficients of the pair correlation function
  // density[k] gives the (un-normalized) coefficient of x^2k in the density
  // correlation[k] gives the (un-normalized) coefficient of x^2k in the correlation function
  // k = 0 to gammaover2 * (N-1)
  
  Density_coefs<T>(int gam_i, long int n_i, unifunction<T>* f_i);  // constructor
  //  Density_coefs<T>(const Density_coefs<T>& a); // copy constructor use default
  Density_coefs<T>(const Density_coefs<T>* dens_ptr); // copy only the functions pointers, but the density is initialized to 0.0
  Density_coefs<T>& operator+= (const Density_coefs<T>& a);
  Density_coefs<T> operator+ (const Density_coefs<T>& a) const;
  void clear();
  void update(const T& prefactor, cpartition& part);
  // prints the coefficients of dens = density or correlation
  // if math_mode == false : output a table
  // if math_mode == true  : output a Mathematica readable form of the polynomial
  void print_density(std::string name, int prec, T& Z);
  void print_correlation(std::string name, int prec, T& Z);

 protected:
  void print(std::string filename, int prec, std::vector<T> dens, T& Z, std::string dens_name, bool math_mode=false);
  void print(std::ostream& os, int prec, std::vector<T> dens, T& Z, std::string dens_name, bool math_mode=false);

};


// Model that computes the checksum, partition function and the coeffcients of the polynomial part 
// of the a set of the density and correlation function
//
template <class T, class Density_Set = Cdensity_set<T> >
  class Model_density_coefs : public Model_density<T, Density_Set>
{
 public:
 Density_coefs<T> density_coefs;
 Model_density_coefs(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, int max_prec_i);
 Model_density_coefs(const Model_density_coefs& mod); // copy all data from mod
 Model_density_coefs(Model_density_coefs * mod); // copy function pointers from *mod, but densities=0.0
 Model_density_coefs<T, Density_Set>& operator+= (const Model_density_coefs<T, Density_Set>& a);
 Model_density_coefs<T, Density_Set> operator+ (const Model_density_coefs<T, Density_Set>& a) const;
 // update results with data from partition part
 T update(cpartition& part);
 void output_report(process_options& op, std::string mp_type, std::string model_name);
 void print_internal_energy(); // outputs the internal energy, implemented only for the sphere_N model
};




#endif
