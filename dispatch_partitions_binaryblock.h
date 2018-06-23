#ifndef __DISPATCH_PARTITIONS_BINARYBLOCK_H_INCLUDED__
#define __DISPATCH_PARTITIONS_BINARYBLOCK_H_INCLUDED__


#include <fstream>
#include "tbb/tbb.h"
#include "partition.h"
#include "cresults.h"
#include "textslice.h"


class partition_dispatcher : public tbb::filter
{
public:
  partition_dispatcher(size_t block_size_i, long int partition_size_i, std::ifstream* readFile_i);
partition_dispatcher(const partition_dispatcher& pd);
~partition_dispatcher();
void fill_buffer() ;
void* operator () (void* );

protected:
  size_t block_size; // = number of partitions in the block
  long int partition_size;
  cpartition part; // dummy partition to get the sizeof(coef) etc...
  size_t raw_block_size; // = size of the block in bytes
  std::ifstream* readFile;
  TextSlice* block_of_data;
  bool work_done;
};

template <class T>
class single_model_worker : public tbb::filter
{
protected:
  long int n; // partition size=number of particles
  Model<T>* model_data;
  // contains pointer to the functions used to compute partition functions

public:
  single_model_worker(long int partition_size_i,
		      Model<T>* model_i);
  void* operator()(void* item);
};

template <class ModelT>
class model_worker : public tbb::filter
{
protected:
  long int n; // partition size=number of particles
  ModelT* model_data;
  // contains pointer to the functions used to compute partition functions

public:
  model_worker(long int partition_size_i,
		      ModelT* model_i);
  void* operator()(void* item);
};

/*
class compute_worker : public tbb::filter
{
protected:
  long int n; // partition size=number of particles
  // pointer to the functions used to compute the scalar long double quantities
  std::vector<unifunction<long double> > *unifnc_ld;
  // pointer to the functions used to compute the scalar mp_real quantities
  std::vector<unifunction<mp_real> > *unifnc_mp;
  // pointer to the functions used to compute the long double densities
  std::vector<function<long double> > *fnc_ld;
  // pointer to the functions used to compute the mp_real densities
  std::vector<function<mp_real> > *fnc_mp;

public:
  compute_worker(long int partition_size_i,
		 std::vector<unifunction<long double> > *unifnc_ld_i,
		 std::vector<unifunction<mp_real> > *unifnc_mp_i,
		 std::vector<function<long double> > *fnc_ld_i,
		 std::vector<function<mp_real> > *fnc_mp_i);
  void* operator()(void* item);
};

class consolidate_results : public tbb::filter
{
protected:
  Cresults* results;
public:
  consolidate_results(Cresults* results_i); 
  void* operator() (void*);
};
*/

template <class T>
class consolidate_single_model : public tbb::filter
{
protected:
  Model<T>* results;
public:
  consolidate_single_model(Model<T>* results_i); 
  void* operator() (void*);
};

template <class ModelT>
class consolidate_model : public tbb::filter
{
protected:
  ModelT* results;
public:
  consolidate_model(ModelT* results_i); 
  void* operator() (void*);
};



#endif
