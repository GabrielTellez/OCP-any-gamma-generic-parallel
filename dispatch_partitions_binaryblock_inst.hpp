#include "arprec/mp_real.h"
#include "mpreal.h"
#include <gmpxx.h>

template class single_model_worker<mp_real>;
template class consolidate_single_model<mp_real>;
template class single_model_worker<mpf_class>;
template class consolidate_single_model<mpf_class>;
template class single_model_worker<mpfr::mpreal>;
template class consolidate_single_model<mpfr::mpreal>;
template class single_model_worker<long double>;
template class consolidate_single_model<long double>;

template class model_worker<Model_unidensity<long double> >;
template class consolidate_model<Model_unidensity<long double> >;
template class model_worker<Model_unidensity<mp_real> >;
template class consolidate_model<Model_unidensity<mp_real> >;

template class model_worker<Model_density<long double> >;
template class consolidate_model<Model_density<long double> >;
template class model_worker<Model_density<mp_real> >;
template class consolidate_model<Model_density<mp_real> >;

template class model_worker<Model_density_coefs<long double> >;
template class consolidate_model<Model_density_coefs<long double> >;
template class model_worker<Model_density_coefs<mp_real> >;
template class consolidate_model<Model_density_coefs<mp_real> >;
