#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <time.h>
#include "tbb/tbb.h"
#include <arprec/mp_real.h>
// #include "mpreal.h"
#include <gmpxx.h>
// #include <boost/ptr_container/ptr_vector.hpp>
// #include <boost/bind.hpp>


bool debug=false;

#include "partition.h"

#include "dispatch_partitions_binaryblock.h"
#include "cresults.h"
#include "process_options.h"
#include "ocp_functions.h"
#include "ocp_models.h"
#include "run_ocp_model.h"

using namespace std;
using namespace tbb;

int main(int argc, char* argv[])
{  


  unsigned int blocksize;

  long int n=0;
  int gammaover2=0;
  bool binmode=true;
  size_t tokens;
  std::string infile;

  process_options op(argc,argv);
  //  option_processer.check_valid();
  tokens=op.tokens;
  blocksize=op.blocksize;
  infile=op.infile;
  debug=op.debug;

  if (debug) {
    std::cout << "Processed options: \n";
    op.print_options();
  }
  // binmode ignored: always use binary mode
  op.binmode=true;
  binmode=true;
  // Initialize the arbitrary precision libraries to opmax_prec digits
  mp::mp_init(op.max_prec+2);
  mp::mpsetprec(op.max_prec);
  mp::mpsetoutputprec(op.max_prec);
  // The GNU mp library expects the precisions in bits rather than digits
  // From arprec doc: "one word contains 48 bits, or about 14.44 digits"
  mpf_set_default_prec(48.0*op.max_prec/14.44);

  if (op.binmode) {
    op.open_readFile();  // read n and gammaover2 from binary file
    op.close_readFile();
    n=op.n;
    gammaover2=op.gammaover2;
  }
  if(debug) std::cout << "Initializing factorials table... ";
  if(!cpartition::init_factorials(n,gammaover2)) {
    std::cerr<<"Cannot initialize factorial tables. Aborting.\n";
    exit(1);
  }
  if(debug)  std::cout << "done.\n";

  std::cout << "Computing OCP partition function and density in parallel using pipeline TBB library.\n";
  std::cout << "N = " << n << endl; ;
  std::cout << "Gamma = " << 2*gammaover2 << endl;
  std::cout << "Reading partitions in " << (binmode?"binary":"text") <<" mode.\n";
  std::cout << "Block size: " << blocksize << "\n";
  std::cout << "Max tokens in flight: " << tokens << "\n";
  if (op.use_arprec) {
    cout << "Requested arbitrary precision: " << op.max_prec << " digits.\n";
  }
  std::cout << "\n";
  if (op.use_arprec) {
    ocp_model<mp_real> ocp(op);
    if(!op.compute_coefs) {
      run_ocp<mp_real> run_calculations(op, "mp_real",&ocp);
    } else {
      run_ocp<mp_real, Model_density_coefs<mp_real> > run_calculations(op, "mp_real",&ocp);
    }
    std::cout << "\n";
  } else {
    op.max_prec=14;
    ocp_model<long double> ocp(op);
    if(!op.compute_coefs) {
      run_ocp<long double> run_calculations(op, "long_double",&ocp);
    } else {
      run_ocp<long double, Model_density_coefs<long double> > run_calculations(op, "long_double",&ocp);
    }
    std::cout << "\n";
  }

  mp::mp_finalize();
  return 0;
}

