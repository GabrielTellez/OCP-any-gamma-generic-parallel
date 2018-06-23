#include <assert.h>
#include "dispatch_partitions_binaryblock.h"
#include "cresults.h"

using namespace tbb;

partition_dispatcher::partition_dispatcher(size_t block_size_i,
					   long int partition_size_i, std::ifstream* readFile_i) 
  : filter(serial_in_order),
    block_size(block_size_i),
    partition_size(partition_size_i),
    part(partition_size),
    raw_block_size(block_size_i*(sizeof(part.coef)+partition_size*sizeof(part.partmembers[0]))),
    readFile(readFile_i),
    block_of_data(TextSlice::allocate(raw_block_size)), work_done(false)
{ }
partition_dispatcher::partition_dispatcher(const partition_dispatcher& pd)
    :     filter(serial_in_order), block_size(pd.block_size),
	  partition_size(pd.partition_size),
	  raw_block_size(pd.raw_block_size),
	  readFile(pd.readFile),
	  block_of_data(pd.block_of_data), work_done(pd.work_done)
{
  assert(!pd.block_of_data); 
}
partition_dispatcher::~partition_dispatcher()
{
  if  (block_of_data) {
    block_of_data->free();
  }
}

void partition_dispatcher::fill_buffer() 
{
  size_t m = block_of_data->avail();
  // debug
  // std::cout << "avail=" <<m<< " , raw_block_size = " << raw_block_size << "\n";
  if (m < raw_block_size ) {
    std::cerr << "Problem: available space in TextSlice is smaller that raw_block_size\n";
    exit(-1);
  }
  readFile->read( block_of_data->end(), raw_block_size);
  size_t n=readFile->gcount();
  block_of_data->set_end(block_of_data->end()+n);
  if(readFile->eof() || n< raw_block_size) {
    //      std::cout << "End of file \n"; // debug
    work_done=true;
  } else {
    work_done=false;
  }
}
		  
void* partition_dispatcher::operator () (void* ) 
{
  if (work_done) {
    return NULL;
  } else {
    fill_buffer();
    TextSlice* send_block=block_of_data;
    block_of_data = TextSlice::allocate(raw_block_size);
    // debug      std::cout << "Send block of size " << send_block->size() << "\n";
    return send_block;
  }
}

/*
compute_worker::compute_worker(long int partition_size_i,
		 std::vector<unifunction<long double> > *unifnc_ld_i,
		 std::vector<unifunction<mp_real> > *unifnc_mp_i,
		 std::vector<function<long double> > *fnc_ld_i,
		 std::vector<function<mp_real> > *fnc_mp_i)
  : filter(parallel), n(partition_size_i),
    unifnc_ld(unifnc_ld_i),
    unifnc_mp(unifnc_mp_i),
    fnc_ld(fnc_ld_i),
    fnc_mp(fnc_mp_i)
{ }

void* compute_worker::operator()(void* item)
{
  // get the textslice with the partitions
  TextSlice* block_of_data=static_cast<TextSlice*>(item);
  // Initialize partial results at 0.0 value
  Cresults partial_res(unifnc_ld->size(), unifnc_mp->size(), fnc_ld, fnc_mp);
  long double factmult=1.0;
  long double coefz=0.0;
  cpartition part(n);
  size_t raw_partition_size=sizeof(part.coef)+part.partmembers.size()*sizeof(part.partmembers[0]);
  char* block_pointer=block_of_data->begin();
  while(block_pointer<block_of_data->end()) {
    part.read_raw(block_pointer);   
    //    std::cout << "Read partition " << part <<" \n" << std::flush; // debug
    factmult=part.mult_factorial_products<long double>();
    coefz=part.coef*part.coef/factmult;
    partial_res.checksum+=coefz;
    for (size_t i=0; i< partial_res.scalar_ld.size(); i++) {
      // compute long double scalars
      partial_res.scalar_ld[i]+=coefz*part.product<long double>(1,n,&(unifnc_ld->at(i)));
    }
    for (size_t i=0; i< partial_res.scalar_ld.size(); i++) {
      partial_res.scalar_mp[i]+=mp_real(coefz)*part.product<mp_real>(1,n,&(unifnc_mp->at(i)));
    }
    block_pointer+=raw_partition_size;
    // TODO: compute densities

  }
  Cresults* res_out=new Cresults(partial_res);

  block_of_data->free();
  return res_out;
  
}

consolidate_results::consolidate_results(Cresults* results_i) 
  : filter(serial_out_of_order), results(results_i)
{}

void* consolidate_results::operator() (void* item)
{
  Cresults* partial=static_cast<Cresults*> (item);
  *results+=*partial;
  delete partial;
  return NULL;
}
*/

template <class T>
single_model_worker<T>::single_model_worker(long int partition_size_i,
					    Model<T>* model_i)
  : filter(parallel), n(partition_size_i), model_data(model_i)
{ }

template <class T>
void* single_model_worker<T>::operator()(void* item)
{
  // get the textslice with the partitions
  TextSlice* block_of_data=static_cast<TextSlice*>(item);
  // Initialize partial results at 0.0 value
  Model<T> partial_res(model_data); // initialize empty results from the template model_data
  long double factmult=1.0;
  long double coefz=0.0;
  cpartition part(n);
  size_t raw_partition_size=sizeof(part.coef)+part.partmembers.size()*sizeof(part.partmembers[0]);
  char* block_pointer=block_of_data->begin();
  while(block_pointer<block_of_data->end()) {
    part.read_raw(block_pointer);   
    //    std::cout << "Read partition " << part <<" \n" << std::flush; // debug
    factmult=part.mult_factorial_products<long double>();
    coefz=part.coef*part.coef/factmult;
    partial_res.update(coefz,part);
    block_pointer+=raw_partition_size;
  }
  Model<T>* res_out=new Model<T>(partial_res);
  block_of_data->free();
  return res_out;
}
					 					 
template <class T>				       
consolidate_single_model<T>::consolidate_single_model(Model<T>* results_i) 
  : filter(serial_out_of_order), results(results_i)
{}

template <class T>
void* consolidate_single_model<T>::operator() (void* item)
{
  Model<T>* partial=static_cast<Model<T>* > (item);
  *results+=*partial;
  delete partial;
  return NULL;
}


/* model_worker class */

template <class ModelT>
model_worker<ModelT>::model_worker(long int partition_size_i,
					    ModelT* model_i)
  : filter(parallel), n(partition_size_i), model_data(model_i)
{ }

template <class ModelT>
void* model_worker<ModelT>::operator()(void* item)
{
  // get the textslice with the partitions
  TextSlice* block_of_data=static_cast<TextSlice*>(item);
  // Initialize partial results at 0.0 value
  ModelT partial_res(model_data); // initialize empty results from the template model_data
  cpartition part(n);
  size_t raw_partition_size=sizeof(part.coef)+part.partmembers.size()*sizeof(part.partmembers[0]);
  char* block_pointer=block_of_data->begin();
  while(block_pointer<block_of_data->end()) {
    part.read_raw(block_pointer);   
    //    std::cout << "Read partition " << part <<" \n" << std::flush; // debug
    partial_res.update(part);
    block_pointer+=raw_partition_size;
  }
  ModelT* res_out=new ModelT(partial_res);
  block_of_data->free();
  return res_out;
}
					 					 
template <class ModelT>				       
consolidate_model<ModelT>::consolidate_model(ModelT* results_i) 
  : filter(serial_out_of_order), results(results_i)
{}

template <class ModelT>
void* consolidate_model<ModelT>::operator() (void* item)
{
  ModelT* partial=static_cast<ModelT* > (item);
  *results+=*partial;
  delete partial;
  return NULL;
}


// Explicit instantiate template classes 
#include "dispatch_partitions_binaryblock_inst.hpp"
