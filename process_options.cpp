#include <iostream>
#include <cstdlib>
#include "process_options.h"
#include "boost/program_options.hpp"
#include "ocp_models.h"  // TOPINDEX_OF_MODELS defined there

process_options::process_options(int argc, char* argv[])
  : debug(false), blocksize(1000), tokens(24), binmode(true), n(0), gammaover2(0), infile(""), use_arprec(false), max_prec(200), model(0), compute_coefs(false)
{
  int bigGamma=0;
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produces this help message")
    ("block-size,b", po::value<unsigned int>(&blocksize)->default_value(1000), "number of partitions in a block that are processed by a thread")
    ("tokens,t", po::value<size_t>(&tokens)->default_value(24), "number of tokens (blocks) in flight")
    ("in-file,i", po::value<std::string>(&infile), "input file with the partitions and coefficients data")
    ("debug,d", "debug mode")
    ("text-mode,x", "input file is in text mode")
    ("Gamma,G", po::value<int>(&bigGamma), "value of the coupling constant Gamma")
    ("particle-number,n", po::value<long int>(&n), "number of particles") 
    ("use-arprec", "use the arbitrary precision library (ARPREC) for calculations")
    ("prec", po::value<int>(&max_prec)->default_value(200), "precision in digits for arbitrary precision calculations") 
    ("model,m", po::value<unsigned int>(&model)->default_value(0), "model to be computed. Use --print-models to see available models.")
    ("p1,p", po::value<long double>(&param1)->default_value(0.0), "extra paramater 1 needed for some models. Use --print-models for details") 
    ("p2", po::value<long double>(&param2)->default_value(0.0), "extra paramater 2 needed for some models. Use --print-models for details") 
    ("print-models", "prints available models")
    ("compute-coefs", "compute the coefficients of the polynomial part of the density and correlation function")
    ("compute-density", po::value<std::vector<std::string> >(&density_domains)->multitoken(),
        "compute the density in the range [xmin, xmax] with size points. Argument arg should be given as xmin,xmax,size " );

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cerr << "Error: " << e.what() << "\n";
    std::cerr << desc << "\n";
    exit(1);
  }
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(0);
  }

  if (vm.count("debug")) {
    debug = true;
  } else {
    debug = false;
  }
  if (vm.count("compute-coefs")) {
    compute_coefs = true;
  } else {
    compute_coefs = false;
  }
  if (vm.count("text-mode")) {
    binmode=false;
  } else {
    binmode=true;
  }
  if (vm.count("use-arprec")) {
    use_arprec=true;
  } else {
    use_arprec=false;
  }

  if (vm.count("Gamma")) {
    if(bigGamma % 2 !=0 ) {
      std::cerr << "Gamma = " << bigGamma << " is not an even integer\nAborting...\n";
      //      std::cerr << desc << "\n";
      exit(1);
    }
    gammaover2=bigGamma/2;
  }
  if ( model > TOPINDEX_OF_MODELS ) {
    std::cout << "Model number "<< model << " out of range.\n";
    ocp_model<long double>::print_models();
    exit(1);
  }
  if ( vm.count("print-models") ) {
    ocp_model<long double>::print_models();
    exit(0);
  }
  if (infile=="" ) {
    std::cerr << "Provide an input filename with -i option\n";
    //    std::cout << desc << "\n";
    exit(1);
  }
}


void process_options::print_options()
{
  std::cout << "debug = " << debug << std::endl;
  std::cout << "blocksize = " << blocksize << std::endl;
  std::cout << "tokens = " << tokens << std::endl;
  std::cout << "binmode = " << binmode << std::endl;
  std::cout << "n = " << n << std::endl;
  std::cout << "gammaover2 = " << gammaover2 << std::endl;
  std::cout << "infile = " << infile << std::endl;
  std::cout << "compute_coefs = " << compute_coefs << std::endl;
  for(std::vector<std::string>::iterator it = density_domains.begin() ;
      it != density_domains.end(); it++) {
    std::cout << "compute density = " << *it << std::endl;
  }
  std::cout << "param1 = " << param1 << std::endl;
}

void process_options::open_readFile()
{
  if (binmode) {
    readFile.open(infile.c_str(), std::ios::in | std::ios::binary); 
    if (readFile.fail()) {
      std::cout << "\nFailed to open binary input file " << infile << "\n";
      exit(1);
    }
    readFile.read(reinterpret_cast<char*>(&n), sizeof(n));
    readFile.read(reinterpret_cast<char*>(&gammaover2), sizeof(gammaover2));
  } else {
    readFile.open(infile.c_str(), std::ios::in ); 
    if (readFile.fail()) {
      std::cout << "\nFailed to open text input file " << infile << "\n";
      exit(1);
    }
  }
}
void process_options::close_readFile()
{
  readFile.close();
}

