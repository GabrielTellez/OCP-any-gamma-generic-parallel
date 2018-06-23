#ifndef __RUN_OCP_MODELS_H_INCLUDED__
#define __RUN_OCP_MODELS_H_INCLUDED__ 

#include <string>
#include <time.h>
#include <boost/bind.hpp>
#include "tbb/tbb.h"
#include "process_options.h"
#include "cresults.h"
using namespace tbb;

template <class T, class ModelT = Model_density<T> >
class run_ocp
{
 public:
  template <class OCPMODEL> run_ocp(process_options& op, std::string mp_type, OCPMODEL* );
  //  Model_density<T> results;
 private:
 run_ocp( const run_ocp<T, ModelT>& other ); // non construction-copyable
 run_ocp& operator=( const run_ocp<T, ModelT>& ); // non copyable
};

template <class T, class ModelT>
template <class OCPMODEL>
  run_ocp<T, ModelT>::run_ocp(process_options& op, std::string mp_type, OCPMODEL* mod)
{
  long int n=0;
  int gammaover2=0;
  op.open_readFile();
  n=op.n;
  gammaover2=op.gammaover2;
  
  if(op.debug) std::cout << "Creating tables of precomputed functions... " << mp_type << " precision..." << std::flush ;
  unifunction<T> unif_T_ocp((n-1)*gammaover2, 
			    boost::bind(&OCPMODEL::f,mod, _1));
  if (op.debug) std::cout << "done\n";
  if (op.debug) std::cout << "Creating initial Model<" << mp_type << ">... " << std::flush ;
  ModelT results(n, gammaover2, &unif_T_ocp, op.max_prec);
  if (op.debug) std::cout << "done\n";
  // decide from options op if it is needed to compute the density or its coefficients
  if (op.debug) std::cout << "Creating tables of functions to compute the density... " << std::flush ;
  function_set<T> f_dens_tables(op.density_domains, 
				boost::bind(&OCPMODEL::f_dens,mod,_1,_2), mod->get_name(), gammaover2, n, mp_type);
  if (op.debug) std::cout << "done\nInitializing densities to zero... " << std::flush ;
  results.initialize_density_set(f_dens_tables, mod->get_comment_line());
  if (op.debug) std::cout << "done\n";
  if (op.debug) std::cout << "Creating <" << mp_type << "> pipeline and filters... " << std::flush ;
  tbb::pipeline pipeline;
  partition_dispatcher pd(op.blocksize,n,&op.readFile);
  pipeline.add_filter(pd);
  model_worker<ModelT> chw(n, &results);
  pipeline.add_filter(chw);
  consolidate_model<ModelT> cons(&results);
  pipeline.add_filter(cons);
  if (op.debug) std::cout << "done.\n";
  if (op.debug) std::cout << "Running the <" << mp_type << "> pipeline... " << std::flush ;
  tick_count t0 = tick_count::now();
  pipeline.run(op.tokens);  
  tick_count t1 = tick_count::now();
  if (op.debug) std::cout << "done.\n" << std::flush;
  op.close_readFile();
  if (op.debug) std::cout << "Normalizing the density... ";
  results.normalize_density(boost::bind(&OCPMODEL::pre_f, mod, _1), 1);
  if (op.debug) std::cout << "done.\n" << std::flush;
  results.free_energy=mod->free_energy(results.Z_mp);
  
  double tiempo=(t1-t0).seconds();
  std::cout << "Computation for model " << mod->get_name() << " <" << mp_type << "> took "<< tiempo <<" seconds.\n";
  clock_t tiempotot=clock();
  std::cout << "Total time "<< tiempotot <<" clicks="<< ((float)tiempotot)/CLOCKS_PER_SEC <<" secs." << std::endl;
  
  results.output_report(op, mp_type, mod->get_name());

}




#endif
