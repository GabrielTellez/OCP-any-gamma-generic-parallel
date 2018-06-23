/* Partition class implementation */


#include "partition.h"
#include <cerrno>
#include <typeinfo>
#include <iomanip>
#include <math.h>

// #include <algorithm>
// #include <exception>


#include "cresults.h"


std::vector<long double> cpartition::factorial; // tables of factorials to compute multiplicity factorial


std::ostream& operator<<(std::ostream& output, cpartition& part)
{
  part.print(output);
  return output;
}


cpartition::cpartition(long int n)
// n= number of elements of the partition
{ 
  if(n>0) {
    partmembers=std::vector<long int>(n,0);
    coef=0;
  } else {
    std::cerr << "The size of the partition," << n << ", must be positive. Creating a partition of size 1.\n";
    partmembers=std::vector<long int>(1,0);
    coef=0;
  }      
}

// read partition and coefficient from a string formated as
// coef [mu1, mu2, ... ,mu_n]
bool cpartition::read(std::string str, bool clear_p=true)
{
  size_t pos1;
  size_t pos2;
  std::string sub("");
  int endloop=0;
  std::vector<long int>::iterator it;
  bool sorted;

  pos1=str.find("[");
  if(pos1==std::string::npos) {
    std::cerr << "Error reading partition, badly formatted string='" << str <<"'\n";
    return false;
  }
  std::string str_coef=str.substr(0,pos1-1);
  if(clear_p) clear();
  errno=0;
  char *endptr;
  coef=strtold(str_coef.c_str(),&endptr);
  if (errno == ERANGE) {
    std::cerr<<"Error reading coeficient of the partition: out of range\n";
    return false;
  }
  pos2=pos1;
  for(it=partmembers.begin() ; it!=partmembers.end() ; it++) {
    *it=0;
  }
  it=partmembers.begin();
  sorted=true;
  do {
    pos2=str.find(",",pos1+1);
    if(pos2==std::string::npos) {
      endloop=1;
      pos2=str.find("]",pos1+1);
    }
    sub=str.substr(pos1+1,pos2-pos1);
    *it=atol(sub.c_str());
    if(it!=partmembers.begin()) {
      if(*(it-1)<*it) {
	sorted=false;
	// The partition is not ordered
      }
    }
    pos1=pos2;
    if(!endloop) {
      if(it!=partmembers.end()) {
	it++;
      } else {
	endloop=true;
	std::cerr << "Partition larger than " << partmembers.size() <<", truncating it.\n";
      }
    }
  } while (!endloop);
  if(!sorted) {
    //    std::cerr << "Ordering partition ... \n";
    sort(); 
  }
  return true;
}
bool cpartition::read(std::istream& is, bool clear_p=true)
{
  std::string str;
  if(getline(is,str)) {
    read(str, clear_p);
    return true;
  } else {
    return false;
  }
}

void cpartition::print(std::ostream& os)
{
  std::vector<long int>::iterator it=partmembers.begin();
  os << std::setiosflags(std::ios::fixed) << std::setprecision(0);
  os << coef << " [";
  while(it!=partmembers.end()) {
    os << *it ;
    it++;
    if(it!=partmembers.end()) {
      os << ",";
    }
  }
  os << "]";
}

// returns the sum of the partition
long double cpartition::sum()
{
  long int suma=0;
  std::vector<long int>::iterator it=partmembers.begin();
  while(it!=partmembers.end()) {
    suma+= *it ;
    it++;
  }
  return suma;
}

// returns the maximum integer in the partition
long int cpartition::max(bool ordered)
{
  long int maxi=0;
  if(!ordered) {
    // Partition not ordered, slow code.
    std::vector<long int>::iterator it=partmembers.begin();
    while(it!=partmembers.end()) {
      if (maxi < *it ) maxi = *it;
      it++;
    }
  } else {
    // ordered partition, the maximum is the first term
    maxi=partmembers[0];
  }
  return maxi;
}

long  int cpartition::multiplicity(long int i, bool ordered)
{
  long int maxi=max(ordered);
  if(i>maxi) {
    return 0;
  } else {
  std::vector<long int>::iterator it=partmembers.begin();
  long int mult=0;
  while(it!=partmembers.end()) {
    if (i== *it ) mult++;
    it++;
  }
  return mult;
  }
}

  bool cpartition::inpartition(long int i, bool ordered)
{
  long int maxi=max(ordered);
  bool ok=false;
  if(i>maxi) {
    return false;
  } else {
  std::vector<long int>::iterator it=partmembers.begin();
  while(it!=partmembers.end()) {
    if (i== *it ) { 
      ok=true;
      break;
    }
    it++;
  }
  return ok;
  }
}

long long int integer_factorial(int i)
{
  if (i<=1) {
    return 1;
  } else {
    return i*integer_factorial(i-1);
  }
}

template <class factorial_type>
factorial_type my_factorial(int i)
{
  factorial_type a_factorial_type_var;
  long long int a_llinteger=0;
  double a_double=0;
  long double a_long_double=0;

  if (typeid(a_factorial_type_var)==typeid(a_llinteger)) {
    return integer_factorial(i);
  }
  if (typeid(a_factorial_type_var)==typeid(a_double)) {
    return tgamma((double)i+1.0);
  }
  if (typeid(a_factorial_type_var)==typeid(a_long_double)) {
    return cpartition::factorial[i];
    //    return tgammal((long double)i+1.0);
  }
  std::cerr << "my_factorial: non defined argument type, returning 1.\n";
  return (factorial_type) 1;
}

template long double my_factorial<long double>(int);

template <class factorial_type>
factorial_type cpartition::mult_factorial_products(bool ordered)
{
  factorial_type prod=(factorial_type) 1;
  if(!ordered) {
    // Partition is not ordered, slow algorithm
    long int i;
    long int maxi=max();
    for(i=1;i<=maxi;i++) {
      prod=prod*my_factorial<factorial_type>(multiplicity(i));
    }
  } else {
    // Partition is ordered
    std::vector<long int>::iterator it=partmembers.begin();
    factorial_type mult=(factorial_type) 1;
    while(it!=partmembers.end()) {
      if(it!=partmembers.begin()) {
	if(*it==*(it-1)) {
	  mult++;
	  prod=prod*mult;
	} else {
	  mult=(factorial_type) 1;
	}
      }
      it++;
    }
  }
  return prod;
}

template long double cpartition::mult_factorial_products<long double>(bool ordered);


void cpartition::sort()
{
  std::sort(partmembers.begin(), partmembers.end(), std::greater<long int>());
}

/*
template <class bignum_type>
bignum_type cpartition::product(int first, int last, bignum_type (*f) (long int))
{
  bignum_type result=(bignum_type) 1;

  if (first<1) {
    std::cerr << "Partition element out of bounds, first <0\n";
    return result;
  }
  if (last>(int)partmembers.size()) {
    std::cerr << "Partition element out of bounds, last > " << partmembers.size()-1 << "\n";
    return result;
  }
  int i;
  for(i=first-1; i<=last-1 ; i++) {
    result=result*f(partmembers[i]);
  }
  return result;
}
*/

template <class bignum_type, class functionoid> 
bignum_type cpartition::product(int first, int last, functionoid* f)
{
  bignum_type result=(bignum_type) 1;

  if (first<1) {
    std::cerr << "cpartition::product: Partition element out of bounds, first < 0\n";
    return result;
  }
  if (last>(int)partmembers.size()) {
    std::cerr << "cpartition::product: Partition element out of bounds, last > " << partmembers.size()-1 << "\n";
    return result;
  }
  int i;
  for(i=first-1; i<=last-1 ; i++) {
    result=result*(*f)(partmembers[i]);
  }
  return result;
}

template <class bignum_type, class functionoid> 
bignum_type cpartition::sum(int first, int last, functionoid* f, size_t x_i)
{
  bignum_type result=bignum_type(0.0);

  if (first<1) {
    std::cerr << "cpartition::sum: Partition element out of bounds, first <0\n";
    return result;
  }
  if (last>partmembers.size()) {
    std::cerr << "cpartition::sum: Partition element out of bounds, last > " << partmembers.size()-1 << "\n";
    return result;
  }
  for(size_t i=first-1; i<=last-1 ; i++) {
    result=result+(*f)(partmembers[i],x_i);
  }
  return result;
}

bool cpartition::write_binary(std::ofstream& os) {
  if (!partmembers.empty()) {
    os.write(reinterpret_cast<char*>(&coef), sizeof(coef));
    os.write(reinterpret_cast<char*>(&partmembers[0]), partmembers.size() * sizeof(partmembers[0]));
    return !(os.fail());
  } else {
    return false;
  }
}

bool cpartition::read_binary(std::ifstream& is, bool clear_part=true) {
  if (clear_part) clear();
  is.read(reinterpret_cast<char*>(&coef), sizeof(coef));
  is.read(reinterpret_cast<char*>(&partmembers[0]), partmembers.size() * sizeof(partmembers[0]));
  return !(is.fail());
}

 bool cpartition::read_raw(char* char_p, bool clear_part) {
  if (clear_part) clear();
  memcpy(&coef,char_p,sizeof(coef));
  memcpy(&partmembers[0], char_p+sizeof(coef), partmembers.size()*sizeof(partmembers[0]));
  return true;
}


class create_partition_exception : public std::exception
{
  virtual const char* what() const throw()
  {
    return "create_partition_exception happened\n";
  }
};


// create a partition from file data
// trows an exception of type create_partition_exception if file could not be read
cpartition::cpartition(std::ifstream& is, // input stream
		       bool binmode, // binary=true or text=false mode
		       long int size) // size of the partition
  : coef(0.0), partmembers(size,0)
{
  bool success=false;
  create_partition_exception e;
  if (binmode) {
    success=read_binary(is, /* clear partition = */ false);
  } else {
    success=read(is, /* clear partition = */ false );
  }
  if(!success) {
    throw e;
    clear();
  }
}

void cpartition::clear() {
  std::vector<long int>::iterator it;
  coef=0;
  for(it=partmembers.begin() ; it!=partmembers.end() ; it++) {
    *it=0;
  }
}


// Instantiate the template to avoid linker errors
#include "partition_inst.hpp"
