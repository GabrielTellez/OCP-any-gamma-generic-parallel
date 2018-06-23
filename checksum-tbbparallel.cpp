#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <time.h>
#include "tbb/tbb.h"
#include <arprec/mp_real.h>
#include "mpreal.h"
#include <gmpxx.h>

bool debug=false;

#include "partition.h"
// #include "factorials_and_gammas_mpreal.cpp"
#include "dispatch_partitions_binaryblock.h"
#include "cresults.h"
#include "process_options.h"
#include "softdisk_functions.h"

using namespace std;
using namespace tbb;

template <class T> void run_calculations(process_options& op, string mp_type="");

template <class T>
void print_report(Model<T>& res, process_options& op) {
  int gammaover2 = op.gammaover2;
  long int n = op.n;
  mp_real chks=gamma(mp_real(gammaover2*n+1))/(pow(gamma(mp_real(gammaover2+1)),mp_real(n))*gamma(mp_real(n+1)));
  cout <<  setiosflags(ios::fixed) << setprecision(0) ;
  //    cout << "Soft disk partition function = " << res.zsoft << "\n";
  cout <<"Long double checksum = "  << res.checksum << endl ;
  cout << "Checksum should be equal to (" << gammaover2 ;
  cout << "*" << n << ")!/( (" << gammaover2 << "!)^" << n <<" *"<< n <<"!) = "  ;
  cout << chks << "\n";

  /*
  if(abs(T(res.checksum)-chks)>=1) {
    cout << "Wrong checksum\n";
  } else {
    cout << "Checksum OK\n";
  }
  */

  cout <<"Partition function (long double precision) = "  << res.Z_ld << endl ;
  cout <<  setiosflags(ios::fixed) << setprecision(0) ;
  cout <<"Partition function (arbitrary precision)   = "  << res.Z_mp << endl ;
}



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
    cout << "Processed options: \n";
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

  {
    ifstream readFile;
    if (op.binmode) {
      readFile.open(op.infile.c_str(), std::ios::in | std::ios::binary); 
      if (readFile.fail()) {
	std::cout << "\nFailed to open input file " << op.infile << "\n";
	exit(1);
      }
      readFile.read(reinterpret_cast<char*>(&n), sizeof(n));
      readFile.read(reinterpret_cast<char*>(&gammaover2), sizeof(gammaover2));
      op.n=n;
      op.gammaover2=gammaover2;
      readFile.close();
    }
  }
  if(debug) std::cout << "Initializing factorials table... ";
  if(!cpartition::init_factorials(n,gammaover2)) {
    cerr<<"Cannot initialize factorial tables. Aborting.\n";
    exit(1);
  }
  if(debug)  std::cout << "done.\n";

  cout << "Computing partition checksum in parallel using pipeline TBB library.\n";
  cout << "N = " << n << endl; ;
  cout << "Gamma = " << 2*gammaover2 << endl;
  cout << "Reading partitions in " << (binmode?"binary":"text") <<" mode.\n";
  cout << "Block size: " << blocksize << "\n";
  cout << "Max tokens in flight: " << tokens << "\n";
  cout << "Requested arbitrary precision: " << op.max_prec << "\n";

  cout << "\n";
  run_calculations<long double>(op, "long double");
  cout << "\n";
  run_calculations<mp_real>(op, "mp_real");
  cout << "\n";
/*
  run_calculations<mpfr::mpreal>(op, "mpfr::mpreal");
  cout << "\n";
  run_calculations<mpf_class>(op, "mpf_class");
  */
  mp::mp_finalize();
  return 0;
}


template <class T> 
void run_calculations(process_options& op, string mp_type)
{
  long int n=0;
  int gammaover2=0;
  ifstream readFile;
  if (op.binmode) {
    readFile.open(op.infile.c_str(), std::ios::in | std::ios::binary); 
    if (readFile.fail()) {
      std::cout << "\nFailed to open input file " << op.infile << "\n";
      exit(1);
    }
    readFile.read(reinterpret_cast<char*>(&n), sizeof(n));
    readFile.read(reinterpret_cast<char*>(&gammaover2), sizeof(gammaover2));
    op.n=n;
    op.gammaover2=gammaover2;
  } else {
    n=op.n;
    gammaover2=op.gammaover2;
  }

  if (op.debug) cout << "Creating tables of precomputed functions: long double precision..." << flush ;
  unifunction<long double> unif_ld_softdisk(n*(n-1)*gammaover2, f_softdisk<long double>);
  if (op.debug) cout << "done.\nCreating tables of precomputed functions: " << mp_type << " precision..." << flush ;
  unifunction<T> unif_mp_softdisk(n*(n-1)*gammaover2, f_softdisk<T>);
  if (op.debug) cout << "done\n";
  if (op.debug) cout << "Creating initial Model<" << mp_type << ">... " << flush ;
  Model<T> softdisk(n, gammaover2, &unif_ld_softdisk, &unif_mp_softdisk);
  if (op.debug) cout << "done\n";

  if (op.debug) cout << "Creating <" << mp_type << "> pipeline and filters... " << flush ;
  tbb::pipeline pipeline;
  partition_dispatcher pd(op.blocksize,n,&readFile);
  pipeline.add_filter(pd);
  single_model_worker<T> chw(n, &softdisk);
  pipeline.add_filter(chw);
  consolidate_single_model<T> cons(&softdisk);
  pipeline.add_filter(cons);
  if (op.debug) cout << "done.\n";

  if (op.debug) cout << "Running the <" << mp_type << "> pipeline... " << flush ;
  tick_count t0 = tick_count::now();
  pipeline.run(op.tokens);  
  tick_count t1 = tick_count::now();
  if (op.debug) cout << "done.\n" << flush;
  readFile.close();
  
  double tiempo=(t1-t0).seconds();
  std::cout << "Computation <" << mp_type << "> took "<< tiempo <<" seconds.\n";
  clock_t tiempotot=clock();
  std::cout << "Total time "<< tiempotot <<" clicks="<< ((float)tiempotot)/CLOCKS_PER_SEC <<" secs." << std::endl;
  
  cout << "Report for <" << mp_type << "> calculations: \n";
  print_report(softdisk, op);



}
