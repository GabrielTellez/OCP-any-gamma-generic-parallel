#ifndef __PROCESS_OPTIONS_H_INCLUDED__
#define __PROCESS_OPTIONS_H_INCLUDED__ 
#include <string>
#include <vector>
#include <fstream>
#include <iostream>


class process_options
{
 public:
  process_options (int argc, char* argv[]);
  bool debug;
  unsigned int blocksize;
  size_t tokens;
  bool binmode;
  long int n;
  int gammaover2;
  std::string infile;
  bool use_arprec;  // use arprec library (true) or long double precision (false)
  int max_prec;    // max precision for arprec library
  unsigned int model; // model to compute = 0 soft disk; 1 hard disk, etc..
  bool compute_coefs; // compute the coefficients of the polynomial part of the density and pair correlation
  std::vector<std::string> density_domains;
  std::ifstream readFile;
  long double param1; // extra parameter for some models
  // W = perimeter (in reduced units) for soft cylinder model
  // Q = colloid charge for Manning 2D
  long double param2; // extra parameter for some models
  // D/R = cell model box size for disk colloid (Manning 2D)
 public:
  void print_options();
  void open_readFile();
  void close_readFile();
};



#endif

