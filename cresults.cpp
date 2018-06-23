#include "cresults.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
// to convert long double to mpf_class
#include "mygmpxx.h"
// Need to compute checksum directly for comparison
#include <arprec/mp_real.h>
#include "ocp_functions.h"

class set_Domain_exception : public std::exception
{
  virtual const char* what() const throw()
  {
    return "set_Domain_exception happened\n";
  }
};


template <class T>
Domain<T>::Domain()
  :  xmin(0.0), xmax(0.0), xsize(2)
{}

template <class T>
Domain<T>::Domain(const T& xmin_i, const T& xmax_i, const size_t  xsize_i) 
  :  xmin(xmin_i), xmax(xmax_i), xsize(xsize_i)
{
  if(size() <2 ) {
    std::cerr << "Domain size has to be >=2, setting it to 2\n";
    xsize=2;
  } 
}
template <class T>
Domain<T>::Domain(const Domain<T> & dom) 
  :  xmin(dom.min()), xmax(dom.max()), xsize(dom.size())
{}
template <class T>
Domain<T>::Domain(const std::string s )
  :  xmin(0.0), xmax(0.0), xsize(2)
{
  set(s);
}

/* Declared inline
template <class T>
T Domain<T>::min() const { return xmin; } 
template <class T>
T Domain<T>::max() const { return xmax; }
template <class T>
size_t Domain<T>::size() const { return xsize; }
template <class T>
T Domain<T>::step() const { return (max()-min())/(size()-1); }
*/
  // Gives the point x corresponding to the index i 
template <class T>
T Domain<T>::x(size_t i) const 
{
  if( i >= size() ) {
    std::cerr << "Domain: warning x out of bounds\n";
  }
  return min()+T(i)*step();
}
template <class T>
bool Domain<T>::operator==(const Domain<T>& dom) const
{
  return ( min()==dom.min() )&&( max()==dom.max() )&&( size()==dom.size());
} 
template <class T>
void Domain<T>::set(std::string str)
{
  //  set_Domain_exception e;
  long double x;
  size_t s;
  size_t pos1=0;
  size_t pos2=str.find(",",pos1);
  if(pos2>=str.length()) {
    std::cerr << "Badly formated density Domain: could not read xmin\n";
    exit(1);
  }
  {
    std::stringstream ss(str.substr(pos1,pos2-pos1));
    ss >> x;
    if(ss.fail()) {
      std::cerr << "Badly formated density Domain: could not read xmin\n";
      exit(1);
    }
  }
  // DEBUG  std::cout << "min=" << x << "\n";
  set_xmin(ld_to_mp<T>(x));
  pos1=pos2+1;
  pos2=str.find(",",pos1);
  if(pos2>=str.length()) {
    std::cerr << "Badly formated density Domain: could not read xmin\n";
    exit(1);
  }
  {
    std::stringstream ss(str.substr(pos1,pos2-pos1));
    ss >> x;
    if(ss.fail()) {
      std::cerr << "Badly formated density Domain: could not read xmax\n";
      exit(1);
    }
  }
  // DEBUG  std::cout << "max=" << x << "\n";
  set_xmax(ld_to_mp<T>(x));
  pos1=pos2+1;
  if(pos1>=str.length()) {
    std::cerr << "Badly formated density Domain: could not read size\n";
    exit(1);
  }
  {
    std::stringstream ss(str.substr(pos1));
    ss >> s;
    if(ss.fail()) {
      std::cerr << "Badly formated density Domain: could not read size\n";
      exit(1);
    }
  }
  // DEBUG  std::cout << "size=" << s << "\n";
  set_size(s);
  //  print_domain(std::cout); DEBUG
}
template <class T>
void Domain<T>::print_domain(std::ostream& os, std::string sep)
{
  os << min() << sep << max() << sep << size() ;
}
template <class T>
std::istream& operator>> (std::istream& in, Domain<T>& dom )
{
  std::string str;
  in >> str;
  dom.set(str);
  return in;
}
template <class T>
std::ostream& operator<< (std::ostream& os, Domain<T>& dom )
{
  dom.print_domain(os);
  return os;
}

/* inline in header
// Doesn't work for T=mp_real
template <class T>
size_t Domain<T>::index(const T& x) const
{
  return size_t((x-this->min())/this->step());
}
*/

template <class T>
DomainId<T>::DomainId(const Domain<T>& dom_i, std::string name_i) 
  : Domain<T> (dom_i), name(name_i)
{}
template <class T>
DomainId<T>::DomainId(const T& xmin_i, const T& xmax_i, const size_t xsize_i, const std::string name_i) 
  : Domain<T> (xmin_i, xmax_i, xsize_i), name(name_i)
{}
template <class T>
DomainId<T>::DomainId(std::string s, std::string name_i) 
  : Domain<T> (s), name(name_i)
{}
template <class T>
DomainId<T>::DomainId (const DomainId<T> & dom_i) 
  : Domain<T> (dom_i), name(dom_i.get_name())
{}
/* inline in header
template <class T>
std::string DomainId<T>::get_name() const
{
  return name;
}
*/
template <class T>
void DomainId<T>::set_name(std::string name_i)
{
  name=name_i;
}

template <class T>
Cdensity_result<T>::Cdensity_result(const DomainId<T>& dom_i, std::string comment_i) 
  : DomainId<T> (dom_i), density(dom_i.size(),0.0), comment(comment_i)
{ }

template <class T>
Cdensity_result<T>::Cdensity_result(const Cdensity_result<T>& a) 
  : DomainId<T> (a), density(a.density), comment(a.comment)
{ }
template <class T>
  T Cdensity_result<T>::operator[](size_t i) { return density[i]; }
template <class T>
bool Cdensity_result<T>::compare_type(const Cdensity_result<T>& a) const
{
  //returns true if a and *this are of the same type: same size same xmin and xmax.
  return *this==a;
}   
template <class T>
Cdensity_result<T> Cdensity_result<T>::operator+ (const Cdensity_result<T>& a) const
{
  if(!compare_type(a)) {
    std::cerr << "Cannot sum densities with different size or different xmin or xmax\n";
    return *this;
  }
  Cdensity_result res(*this);
  for(size_t i=0;i < this->size();i++) {
    res.density[i]=res.density[i]+a.density[i];
  }
  return res;   //returns *this + a,  *this is not modified
}
template <class T>
Cdensity_result<T>& Cdensity_result<T>::operator+= (const Cdensity_result<T>& a)
{
  if(!compare_type(a)) {
    std::cerr << "Cannot sum densities with different size or different xmin or xmax\n";
    return *this;
}
  for(size_t i=0;i< this->size();i++) {
    density.at(i)+=a.density[i];
  }
  return *this;
}
template <class T>
Cdensity_result<T> Cdensity_result<T>::operator/(const T& z) const
{
  Cdensity_result res(*this);
  for(size_t i=0;i< this->size();i++) {
    res.density.at(i)=density.at(i)/z;
  }
  return res;
}
template <class T>
Cdensity_result<T> Cdensity_result<T>::normalize(const T& z)
{
  for(size_t i=0;i<this->size();i++) {
    density.at(i)=density.at(i)/z;
  }
  return *this;
}
template <class T>
void Cdensity_result<T>::clear()
{
  for(size_t i=0;i<this->size();i++) {
    density.at(i)=0.0;
  }
}

template <class T>
void Cdensity_result<T>::print(std::ostream& os, int prec)
{
  if (comment!="") {
    os << "# " << comment << "\n";
  }
  for(size_t i=0 ; i< this->size(); i++) {
    os <<  setiosflags(std::ios::scientific) << std::setprecision(prec) ;
    os << this->x(i) << "     " << density.at(i) << "\n";
  }  
}
template <class T>
void Cdensity_result<T>::print(std::string filename, int prec)
{
  filename+= "_";
  filename+=this->get_name();
  filename+=".dat";
  std::fstream fstr(filename.c_str(), std::fstream::out); 
  print(fstr, prec);
  fstr.close();
}
				     
template <class T>
void Cdensity_result<T>::update(const T& prefactor, cpartition& part, function<T>* fnc_dens)
{
  int n=part.partmembers.size();
  for( size_t x_index = 0 ; x_index < density.size() ; x_index++ ) {
    density.at(x_index)+= prefactor * part.sum<T>(1,n,fnc_dens,x_index);
  }
}
template <class T>
Cdensity_result<T> Cdensity_result<T>::operator*=(T (*fnc)(const T&))
{
  for (size_t x_index = 0 ; x_index < density.size() ; x_index++ ) {
    (density.at(x_index))*=fnc(this->x(x_index));
  }
  return *this;
}

// Cdensity class (BEGIN)
template <class  T>
Cdensity<T>::Cdensity(const Cdensity& a)
  : Cdensity_result<T> (a), fnc_dens(a.fnc_dens)
{ }
template <class  T>
Cdensity<T>::Cdensity(const Cdensity* dens_ptr)
  : Cdensity_result<T> (static_cast<DomainId<T> >(*dens_ptr)), fnc_dens(dens_ptr->fnc_dens)
{ }
template <class T>
void Cdensity<T>::update(const T& prefactor, cpartition& part)
{
  Cdensity_result<T>::update(prefactor, part, fnc_dens);
}

// Cdensity class (END)

/*
Cresults::Cresults (size_t num_scalar_ld, size_t num_scalar_mp) 
  : checksum(0.0), scalar_ld(num_scalar_ld,0.0), scalar_mp(num_scalar_mp,0.0) 
{}
// Initialize Cresults with DomainId
Cresults::Cresults (size_t num_scalar_ld, size_t num_scalar_mp, 
		    std::vector<DomainId<long double> > density_ld_dom, 
		    std::vector<DomainId<mp_real> > density_mp_dom)
  : checksum(0.0), scalar_ld(num_scalar_ld,0.0), scalar_mp(num_scalar_mp,0.0)
{
  for (size_t i=0; i< density_ld_dom.size() ; i++) {
    density_ld.push_back(Cdensity_result<long double>(density_ld_dom[i]));
  }
  for (size_t i=0; i< density_mp_dom.size(); i++ ) {
    density_mp.push_back(Cdensity_result<mp_real>(density_mp_dom[i]));
  }
}
// Initialize Cresults with function pointer (using their DomainId)
Cresults::Cresults (size_t num_scalar_ld, size_t num_scalar_mp, 
		    std::vector<function<long double> >* func_ld, 
		    std::vector<function<mp_real> >* func_mp)
  : checksum(0.0), scalar_ld(num_scalar_ld,0.0), scalar_mp(num_scalar_mp,0.0)
{
  for (size_t i=0; i< func_ld->size() ; i++) {
    density_ld.push_back(Cdensity_result<long double>((func_ld->at(i)).dom));
  }
  for (size_t i=0; i< func_mp->size() ; i++) {
    density_mp.push_back(Cdensity_result<mp_real>((func_mp->at(i)).dom));
  }
}


Cresults::Cresults (const Cresults& res) 
  : checksum(res.checksum), scalar_ld(res.scalar_ld), scalar_mp(res.scalar_mp), density_ld(res.density_ld),
    density_mp(res.density_mp)
{ }

bool Cresults::compare_type(const Cresults& a) const 
{
  bool compare =  (scalar_ld.size()==a.scalar_ld.size() ) && (scalar_mp.size()==a.scalar_mp.size() ) && (density_ld.size() == a.density_ld.size()) && (density_mp.size() == a.density_mp.size()) ;
  if (!compare) {
    return false;
  } else {
    for (size_t i=0; i<density_ld.size(); i++) {
      compare = compare && (density_ld.at(i).compare_type(a.density_ld.at(i)));
      if (!compare) return false;
    }
    for (size_t i=0; i<density_mp.size(); i++) {
      compare = compare && (density_mp.at(i).compare_type(a.density_mp.at(i)));
      if (!compare) return false;
    }
    return compare;
  }
  return compare;
}
*/

// Cdensity_set class (BEGIN)

template <class T>
Cdensity_set<T>::Cdensity_set()
  : std::vector<Cdensity<T> >()
{}
template <class T>
Cdensity_set<T>::Cdensity_set(const Cdensity_set<T>& a)
  : std::vector<Cdensity<T> >(a)
{}
template <class T>
Cdensity_set<T>::Cdensity_set(Cdensity_set<T>* set_ptr)
  : std::vector<Cdensity<T> >(*set_ptr)
{
  // clear the densities
  for( size_t i=0 ; i < set_ptr->size();  i++) {
    this->at(i).clear();
  }
}


template <class T>
Cdensity_set<T> Cdensity_set<T>::operator+ (const Cdensity_set<T>& a) const 
{
  Cdensity_set res(*this);
  res+=a;
  return res;
}

template <class T>
Cdensity_set<T>& Cdensity_set<T>::operator+= (const Cdensity_set<T>& a)
{
  if ( this->size()!=a.size() ) {
    std::cerr << "Cannot sum Cdensity_Sets of different size\n";
    exit(1);
  }
  for (size_t i=0; i<this->size(); i++) {
    (this->at(i))+=a.at(i);
  }
  return *this;
}

template <class T>
void Cdensity_set<T>::print(std::string filename, int prec)
{
  for (size_t i=0; i<this->size(); i++) {
    (this->at(i)).print(filename, prec);
  }
}
template <class T>
void Cdensity_set<T>::update(const T & prefactor, cpartition& part)
{
  for (size_t i=0; i<this->size(); i++) {
    (this->at(i)).update(prefactor, part);
  }
}
template <class T>
void Cdensity_set<T>::initialize_new_density(function<T>* f_i, std::string com)
{
  this->push_back(Cdensity<T>(f_i->dom, f_i, com));
}
template <class T>
void Cdensity_set<T>::initialize_density_set(function_set<T>& ftables, std::string com)
{
  for (size_t i=0; i < ftables.size() ; i++ ) {
    initialize_new_density(ftables[i],com);
  }
}


// Cdensity_set class (END)

//template <class T>
//function<T>::function(long int n_i,DomainId<T> dom_i, const T (*fnc)(unsigned int,const T&))
//  : n(n_i), dom(dom_i), table(n+1, std::vector<T>(dom.size(),0.0))
/*
template <class T, class F>
function<T>::function(long int n_i,DomainId<T> dom_i, F fnc)
  : n(n_i), dom(dom_i), table(n+1, std::vector<T>(dom.size(),0.0))
{
  // initialize the table of precomputed values
  for(unsigned int mu=0; mu<=size_t(n); mu++) {
    for(unsigned int i=0; i< dom.size(); i++) {
      table[mu][i]=fnc(mu, dom.x(i));
    }
  }
}
*/
/* dom.index not defined for T=mp_real
template <class T>
inline
T function<T>::operator() (long int mu, T x) const
{
  // Use the precomputed table to give the value (faster than computing *f directly)
  return table[size_t(mu)][dom.index(x)];
}
*/

template <class T>
inline
T function<T>::operator() (long int mu, size_t index_x_i) const
{
  // Use the precomputed table to give the value (faster than computing *f directly)
  return table[size_t(mu)][index_x_i];
}

template <class T>
inline
T function<T>::operator() (long int mu) const
{
  // Use the precomputed table to give the value (faster than computing *f directly)
  return table[size_t(mu)][index_x];
}

template <class T>
inline
T unifunction<T>::operator() (long int mu) const
{
  // Use the precomputed table to give the value (faster than computing *f directly)
  return table.at(size_t(mu));
}

template <class T> 
unifunction<T>::unifunction(const unifunction<T>& a)
  : table(a.table)
{}


template <class T>
Model<T>::Model(long int n_i, int gam_i, unifunction<long double>* fnc_ld_i, unifunction<T>* fnc_mp_i)
  : n(n_i), gammaover2(gam_i), checksum(0.0), Z_ld(0.0), Z_mp(0.0),
    fnc_ld(fnc_ld_i), fnc_mp(fnc_mp_i)
{}

// copy model data from mod but with empty results = 0.0
template <class T>
Model<T>::Model(Model<T> * mod)
  : n(mod->n), gammaover2(mod->gammaover2), checksum(0.0), Z_ld(0.0), Z_mp(0.0),
    fnc_ld(mod->fnc_ld), fnc_mp(mod->fnc_mp)
{}

template <class T>
Model<T>& Model<T>::operator+= (const Model<T>& a)
{
  checksum+=a.checksum;
  Z_ld+=a.Z_ld;
  Z_mp+=a.Z_mp;
  // TO DO: densities and correlations
  return *this;
}

template <class T>
Model<T> Model<T>::update(long double coefz, cpartition& part)
{
  checksum+=coefz;
  Z_ld+=coefz*part.product<long double>(1,n,fnc_ld);
  Z_mp+=ld_to_mp<T>(coefz)*part.product<T>(1,n,fnc_mp);
  // TO DO: densities and correlations
  return *this;
}


template <class T>
Model<T> Model<T>::operator+ (const Model<T>& a) const
{
  Model res(*this);
  res+=a;
  return res;
}


template <class T>
ModelBase<T>::ModelBase(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, int max_prec_i)
  : n(n_i), gammaover2(gam_i), count(0.0), checksum(0.0), Z_mp(0.0),
    fnc_mp(fnc_mp_i), max_prec(max_prec_i)
{}

// copy all model data from mod 
template <class T>
ModelBase<T>::ModelBase(const ModelBase<T>& mod)
  : n(mod.n), gammaover2(mod.gammaover2), count(mod.count), checksum(mod.checksum), 
    Z_mp(mod.Z_mp), fnc_mp(mod.fnc_mp), max_prec(mod.max_prec)
{}

// copy model data from mod but with empty results = 0.0
template <class T>
ModelBase<T>::ModelBase(ModelBase<T> * mod)
  : n(mod->n), gammaover2(mod->gammaover2), count(0.0), checksum(0.0), Z_mp(0.0),
    fnc_mp(mod->fnc_mp), max_prec(mod->max_prec)
{}

template <class T>
ModelBase<T>& ModelBase<T>::operator+= (const ModelBase<T>& a)
{
  count+=a.count;
  checksum+=a.checksum;
  Z_mp+=a.Z_mp;
  return *this;
}

template <class T>
T ModelBase<T>::update(cpartition& part)
{
  T factmult=ld_to_mp<T>(part.mult_factorial_products<long double>());
  T partcoef = ld_to_mp<T>(part.coef);
  T coefz=partcoef*partcoef/factmult;
  T prefactor=coefz*part.product<T>(1,n,fnc_mp);
  count+=1;
  checksum+=coefz;
  Z_mp+=prefactor;
  return prefactor;
}


template <class T>
ModelBase<T> ModelBase<T>::operator+ (const ModelBase<T>& a) const
{
  ModelBase<T> res(*this);
  res+=a;
  return res;
}

template <class T> 
void ModelBase<T>::print_report(std::ostream& os, process_options& op, std::string mp_type, std::string model_name)
{
  mp_real chks=gamma(mp_real(gammaover2*n+1))/(pow(gamma(mp_real(gammaover2+1)),mp_real(n))*gamma(mp_real(n+1)));
  os << "Report for " <<  model_name << " <"  << mp_type << "> calculations: \n";
  os << "OCP with N = " << n << " particles at Gamma = "<< 2*gammaover2 <<"\n"; 
  os << "Parameter p1 = " << op.param1 << "\n";
  os << "Parameter p2 = " << op.param2 << "\n";
  os <<  setiosflags(std::ios::fixed) << std::setprecision(0) ;
  os << "Total number of partitions = " << count << "\n";
  os <<  setiosflags(std::ios::fixed) << std::setprecision(0) ;
  os << "Checksum = "  << checksum << "\n";
  os << "Checksum should be equal to (" << gammaover2 ;
  os << "*" << n << ")!/( (" << gammaover2 << "!)^" << n <<" *"<< n <<"!) = "  ;
  os << chks << "\n";
  /*
  if(abs(mp_real(res.checksum)-chks)>=1) {
    cout << "Wrong checksum\n";
  } else {
    cout << "Checksum OK\n";
  }
  */
  os <<  setiosflags(std::ios::fixed) << std::setprecision(0) ;
  os << "Partition function (integer part) ";
  if( mp_type != "") {
    os <<"(" << mp_type << " precision) = ";
  } else {
    os << " = ";
  }
  os << Z_mp << "\n";
  os <<  setiosflags(std::ios::scientific) << std::setprecision(max_prec) ;
  os << "Partition function (precision: " << max_prec << ")";
  os << " = ";
  os << Z_mp << "\n";
  os << "Free energy: beta F = "<< free_energy << " + N (1-(Gamma/4))log (rho_b L^2) \n";
}

// class Model_unidensity (BEGIN)
template <class T>
Model_unidensity<T>::Model_unidensity(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, function<T>* fnc_dens_i, int max_prec_i)
  : ModelBase<T>(n_i, gam_i, fnc_mp_i, max_prec_i), density(fnc_dens_i->dom), fnc_dens(fnc_dens_i)
{}

// copy model data from mod but with empty results = 0.0
template <class T>
Model_unidensity<T>::Model_unidensity(Model_unidensity<T> * mod)
  : ModelBase<T>(mod), density(mod->fnc_dens->dom), fnc_dens(mod->fnc_dens) 
{}

template <class T>
Model_unidensity<T>& Model_unidensity<T>::operator+= (const Model_unidensity<T>& a)
{
  ModelBase<T>::operator+= (a);
  density+=a.density;
  return *this;
}

template <class T>
T Model_unidensity<T>::update(cpartition& part)
{
  T prefactor=ModelBase<T>::update(part);
  density.update(prefactor,part,fnc_dens);
  return prefactor;
    // return *this; OLD
}

template <class T>
Model_unidensity<T> Model_unidensity<T>::operator+ (const Model_unidensity<T>& a) const
{
  Model_unidensity<T> res(*this);
  res+=a;
  return res;
}
template <class T>
void Model_unidensity<T>::normalize_density(T (*fnc)(const T&), const T& factor)
{
  density.normalize(fnc, (this->Z_mp));
}
// class Model_unidensity (END)

// class Model_density (BEGIN)
template <class T, class Density_Set>
Model_density<T, Density_Set>::Model_density(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, int max_prec_i)
  : ModelBase<T>(n_i, gam_i, fnc_mp_i, max_prec_i)
{}

// copy all data from mod
template <class T, class Density_Set>
Model_density<T, Density_Set>::Model_density(const Model_density<T, Density_Set>& mod)
  : ModelBase<T>(mod), density(mod.density)
{}
// copy model data from mod but with empty results and densities = 0.0
template <class T, class Density_Set>
Model_density<T, Density_Set>::Model_density(Model_density<T, Density_Set> * mod)
  : ModelBase<T>(mod), density(&(mod->density))
{}

template <class T, class Density_Set>
Model_density<T, Density_Set>& Model_density<T, Density_Set>::operator+= (const Model_density<T, Density_Set>& a)
{
  ModelBase<T>::operator+= (a);
  density+=a.density;
  return *this;
}

template <class T, class Density_Set>
T Model_density<T, Density_Set>::update(cpartition& part)
{
  T prefactor=ModelBase<T>::update(part);
  density.update(prefactor,part);
  return prefactor;
}

template <class T, class Density_Set>
Model_density<T, Density_Set> Model_density<T, Density_Set>::operator+ (const Model_density<T, Density_Set>& a) const
{
  Model_density<T, Density_Set> res(*this);
  res+=a;
  return res;
}

template <class T, class Density_Set>
void Model_density<T, Density_Set>::initialize_new_density(function<T>* f_i, std::string com)
{
  density.initialize_new_density(f_i, com);
}

template <class T, class Density_Set>
void Model_density<T, Density_Set>::initialize_density_set(function_set<T>& f_ptr, std::string com)
{
  density.initialize_density_set(f_ptr, com);
}

template <class T, class Density_Set>
void Model_density<T, Density_Set>::output_report(process_options& op, std::string mp_type, std::string model_name)
{
  ModelBase<T>::print_report(std::cout, op, mp_type,  model_name);
  density.print("density",op.max_prec);
}


// class Model_density (END)

// class Density_coefs

template <class T>
Density_coefs<T>::Density_coefs(int gam_i, long int n_i, unifunction<T>* f_i)
  : gammaover2(gam_i), n(n_i), fnc(f_i), density(gammaover2*(n-1)+1,0.0),
    correlation(gammaover2*(n-1)+1,0.0)
{

}

// copy only the functions pointers, but the density is initialized to 0.0
template <class T>
Density_coefs<T>::Density_coefs(const Density_coefs<T>* dens_ptr)
  : gammaover2(dens_ptr->gammaover2), n(dens_ptr->n), fnc(dens_ptr->fnc), density(gammaover2*(n-1)+1,0.0),
    correlation(gammaover2*(n-1)+1,0.0)
{

}

template <class T>
Density_coefs<T>& Density_coefs<T>::operator+= (const Density_coefs<T>& a)
{
  if ( (gammaover2!=a.gammaover2) || ( n!=a.n ) ) {
    std::cerr << "Cannot sum Density_coefs of different n or gammaover 2\n";
    exit(1);
  }

  for(size_t i=0;i<=size_t(gammaover2*(n-1));i++) {
    density.at(i)+=a.density.at(i);
    correlation.at(i)+=a.correlation.at(i);
  }
  return *this;
}

template <class T>
Density_coefs<T> Density_coefs<T>::operator+ (const Density_coefs<T>& a) const
{
  Density_coefs<T> res(*this);
  res+=a;
  return res;
}

template <class T>
void Density_coefs<T>::clear()
{
  for(size_t i=0;i<=size_t(gammaover2*(n-1));i++) {
    density.at(i)=0.0;
    correlation.at(i)=0.0;
  }
}

template <class T>
void Density_coefs<T>:: update(const T& prefactor, cpartition& part)
{
  if (part.partmembers.size()!=n) {
    std::cerr << "Error: the partition has a size " << part.partmembers.size() << " different from N =" << n << "\n";
    exit(1);
  }
  for(size_t k=0;k < n; k++) {
    long int mu = part.partmembers[k];
    if ( (mu < 0) || (mu>gammaover2*(n-1)) ) {
      std::cerr << "Error: an element of the partition mu(" << k << ")=" << mu << " is out of bounds [0,"<< gammaover2*(n-1) << "]\n";
      exit(1);
    }
    density.at(mu)+=prefactor/(*fnc)(mu);
    if ( (part.partmembers[n-1] == 0) && ( mu != 0) ) {
      correlation.at(mu)+=prefactor/((*fnc)(mu) * (*fnc)(0));
    }
  }
}
template <class T>
void Density_coefs<T>::print(std::ostream& os, int prec, std::vector<T> dens, T& Z, std::string dens_name, bool math_mode)
{
  if (math_mode) {
      os <<  setiosflags(std::ios::fixed) << std::setprecision(prec) ;
      os << dens_name << "u" << n << "[xd_]:=" ;
      for(size_t i=0 ; i< dens.size(); i++) {
	os << dens.at(i) << " * xd^" << i ;
	if (i!= dens.size()-1) {
	  os << " + " ;
	}
      }  
      os << std::endl;
      /*
      os << dens_name << n << "[xd_]:=" << dens_name << "u" << n << "[xd]/" << Z;
      os << std::endl;
      */
      os << dens_name << n << "[xd_]:=" ;
      for(size_t i=0 ; i< dens.size(); i++) {
	os << dens.at(i)/Z << " * xd^" << i ;
	if (i!= dens.size()-1) {
	  os << " + " ;
	}
      }  
  } else {
    os << "#  " << dens_name ;
    os << "# power       unnormalized coefficient     normalized coefficient\n";
    for(size_t i=0 ; i< dens.size(); i++) {
      os <<  setiosflags(std::ios::scientific) << std::setprecision(prec) ;
      os << i << "     " << dens.at(i) << "       " << dens.at(i)/Z <<"\n";
    }  
  }
}
template <class T>
void Density_coefs<T>::print(std::string filename, int prec, std::vector<T> dens, T& Z, std::string name, bool math_mode)
{
  if (math_mode) {
    filename+=".m";
  } else {
    filename+=".dat";
  }
  std::fstream fstr(filename.c_str(), std::fstream::out); 
  print(fstr, prec, dens, Z, name, math_mode);
  fstr.close();
}

template <class T>
void Density_coefs<T>::print_density(std::string name, int prec, T& Z)
{
  std::stringstream filename;
  filename << "density_";
  filename << name << "_G" << 2*gammaover2;
  filename << "_N" << n;
  print(filename.str(), prec, density, Z, "density", false);
  print(filename.str(), prec, density, Z, "density", true);
}

template <class T>
void Density_coefs<T>::print_correlation(std::string name, int prec, T& Z)
{
  std::stringstream filename;
  filename << "correlation_";
  filename << name << "_G" << 2*gammaover2;
  filename << "_N" << n;
  print(filename.str(), prec, correlation, Z, "correlation", false);
  print(filename.str(), prec, correlation, Z, "correlation", true);
}

// class Density_coefs (END)



// class Model_density_coefs (BEGIN)
template <class T, class Density_Set>
Model_density_coefs<T, Density_Set>::Model_density_coefs(long int n_i, int gam_i, unifunction<T>* fnc_mp_i, int max_prec_i)
  : Model_density<T, Density_Set>(n_i, gam_i, fnc_mp_i, max_prec_i), density_coefs(gam_i, n_i, fnc_mp_i)
{}

// copy all data from mod
template <class T, class Density_Set>
Model_density_coefs<T, Density_Set>::Model_density_coefs(const Model_density_coefs<T, Density_Set>& mod)
  : Model_density<T, Density_Set>(mod), density_coefs(mod.density_coefs)
{}
// copy model data from mod but with empty results and densities = 0.0
template <class T, class Density_Set>
Model_density_coefs<T, Density_Set>::Model_density_coefs(Model_density_coefs<T, Density_Set> * mod)
  : Model_density<T, Density_Set>(mod), density_coefs(&(mod->density_coefs))
{}

template <class T, class Density_Set>
Model_density_coefs<T, Density_Set>& Model_density_coefs<T, Density_Set>::operator+= (const Model_density_coefs<T, Density_Set>& a)
{
  Model_density<T, Density_Set>::operator+= (a);
  density_coefs+=a.density_coefs;
  return *this;
}

template <class T, class Density_Set>
T Model_density_coefs<T, Density_Set>::update(cpartition& part)
{
  T prefactor=Model_density<T, Density_Set>::update(part);
  density_coefs.update(prefactor,part);
  return prefactor;
}

template <class T, class Density_Set>
Model_density_coefs<T, Density_Set> Model_density_coefs<T, Density_Set>::operator+ (const Model_density_coefs<T, Density_Set>& a) const
{
  Model_density_coefs<T, Density_Set> res(*this);
  res+=a;
  return res;
}

template <class T, class Density_Set>
void Model_density_coefs<T, Density_Set>::output_report(process_options& op, std::string mp_type, std::string model_name) 
{
  Model_density<T, Density_Set>::output_report(op, mp_type,  model_name);
  // output the coefficients 
  std::string name(model_name);
  name+="_";
  name+=mp_type;
  density_coefs.print_density(name,op.max_prec, this->Z_mp);
  density_coefs.print_correlation(name,op.max_prec, this->Z_mp);
  if(model_name == "sphere_N") {
    print_internal_energy(); // only applies for the sphere_N model 3 so far
  }
}

template <class T, class Density_Set>
void Model_density_coefs<T, Density_Set>::print_internal_energy()
{
  // This is only valid for the sphere N model !!!
 
  T Upp=0.0; // particle-particle energy
  T prefactor=factorial<T>( ((this->n)-1)*(this->gammaover2)+1);
  for(size_t i=1 ; i < density_coefs.correlation.size(); i++) {
    Upp = Upp + density_coefs.correlation.at(i) *  factorial<T>(i)*factorial<T>(((this->n-1)*this->gammaover2-i))*(HarmonicNumber<T>((unsigned int) i)-HarmonicNumber<T>((this->n-1)*this->gammaover2+1)) ;
  }
  Upp=-Upp/(4.0*this->Z_mp* prefactor);
  T U;
  U=Upp-T(this->n*this->n)/4.0;
  U+=T(this->n)*log(T(this->n))/4.0;
  std::cout << "\n";
  std::cout << "Internal energy U/q^2 = " <<  U << " - (1/4) N*log(pi rho_b L^2) \n";
  std::cout << "Internal energy U/q^2 = " <<  U - T(this->n)*log(T(M_PI))/4.0 << " - (1/4) N*log(rho_b L^2) \n";
  std::cout << "Particle-particle energy Upp/q^2 = " << Upp << "\n";
}

// class Model_density_coefs (END)




// explicit instantiate template classes
#include "cresults_inst.hpp"



