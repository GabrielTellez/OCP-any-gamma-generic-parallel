#include "arprec/mp_real.h"
#include "mpreal.h"
#include <gmpxx.h>

// instantiate the (uni)functions classes with T=long double and mp_real
/*
template class function<long double>;
template class function<mp_real>;
*/
template class unifunction<long double>;
template class unifunction<mp_real>;
template class unifunction<mpfr::mpreal>;
template class unifunction<mpf_class>;

template class function<long double>;
template class function<mp_real>;


// instantiate the Model with T=mp_real

template class Model<long double>;
template class Model<mp_real>;
template class Model<mpfr::mpreal>;
template class Model<mpf_class>;

template class ModelBase<long double>;
template class ModelBase<mp_real>;

template class Model_unidensity<long double>;
template class Model_unidensity<mp_real>;

template class Model_density<long double>;
template class Model_density<mp_real>;

// instantiate Domain

template class Domain<long double>;
template class Domain<mp_real>;
template std::istream& operator>> (std::istream&, Domain<long double>&);
template std::istream& operator>> (std::istream&, Domain<mp_real>&);
template std::ostream& operator<< (std::ostream&, Domain<long double>&);
template std::ostream& operator<< (std::ostream&, Domain<mp_real>&);


// instantiate DomainId

template class DomainId<long double>;
template class DomainId<mp_real>;

// instantiate Cdensity_results

template class Cdensity_result<long double>;
template class Cdensity_result<mp_real>;
template class Cdensity<long double>;
template class Cdensity<mp_real>;
template class Cdensity_set<long double>;
template class Cdensity_set<mp_real>;


// instantiate Density_coefs

template class Density_coefs<long double>;
template class Density_coefs<mp_real>;


// instantiate Model_density_coefs

template class Model_density_coefs<long double>;
template class Model_density_coefs<mp_real>;
