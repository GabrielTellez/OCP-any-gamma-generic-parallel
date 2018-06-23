/* Partition class header */
#ifndef __PARTITION_H_INCLUDED__  
#define __PARTITION_H_INCLUDED__   //



#include <string>
#include <vector>
#include <fstream>
#include <iostream>


class cpartition {
 public:
  long double coef; // coefficient in the expansion of the vandermonde
  std::vector<long int> partmembers; // the partition itself
  cpartition(long int n); // create partition with n zeros
 cpartition() 
   : coef(0), partmembers(1,0) {}; // create partition 0^1 with coef=0.
  cpartition(std::ifstream& is, bool binmode, long int size); // create partition of size size and read it from stream is in binary mode if binmode=true, text mode otherwise
  void clear(); // sets coef=0 and partition=0^n

  // read partition and coefficient from a string or stream formatted as
  // coef [mu1, mu2, ... ,mu_n]
  bool read(std::string str,bool);
  bool read(std::istream& is,bool);

  // read partition and coefficient from a ifstream in binary mode
  bool read_binary(std::ifstream& is,bool);
  // read partition and coefficient from a raw block of chars, pointed by char_p: format long double coef, n* long int partition
  bool read_raw(char* char_p, bool clear_part=true);

  //  partition and coefficient to a ofstream in binary mode
  bool write_binary(std::ofstream& os);

  // prints the partition in text readable format
  // coef [mu1, mu2, ... ,mu_n]
  void print(std::ostream& os);

  // returns the sum of the partition
  long double sum();

  // returns the maximum integer in the partition
  long int max(bool ordered=true);

  // returns the multiplicity of integer i in the partition
  long int multiplicity(long int i, bool ordered=true);
  
  // returns the product of the factorials of the multiplicities
  template <class factorial_type> 
    factorial_type  mult_factorial_products(bool ordered=true);

  // tells of the integer i is in the partition
  bool inpartition(long int i, bool ordered=true);

  // orders the partition from largest to smallest;
  void sort();

  // computes the product or sum (=operation) of a function of each
  // element of the partition:

  /* Should work with the generalized functionoid_ptr
  template <class bignum_type> 
    bignum_type product(int first, int last, bignum_type (*f)(long int));
  */
  // computes product of functions of mu_k
  // prod_{k=first}^{last} f(mu_k)
  template <class bignum_type, class functionoid> 
    bignum_type product(int first, int last, functionoid* f);

  // computes sum of functions
  // sum_{k=first}^{last} f(mu_k,x)
  /*
  template <class bignum_type, class functionoid* f> 
    bignum_type sum(int first, int last, functionoid* f);
  */
  template <class bignum_type, class functionoid> 
    bignum_type sum(int first, int last, functionoid* f, size_t x_i);
 public:
  static std::vector<long double> factorial; // tables of factorials to compute multiplicity factorial
  static bool init_factorials(int n, int gammaover2)
  // compute tables of factorials 
  // up to gammaover2* (n-1) + 4
  {
    if ((n<=1)||(gammaover2<1)) {
      std::cerr << "Cannot initialize factorials table: n=" << n << ", Gamma/2=" << gammaover2 <<"\n"; 
      return false;
    }
    factorial=std::vector<long double> (gammaover2*(n-1)+4,1.0);
    for(size_t i=2; i< factorial.size(); i++) {
      factorial[i]=((long double) i)*factorial[i-1];
    }
    return true;
  }
};

std::ostream& operator<<(std::ostream& output, const cpartition& part);


#endif

