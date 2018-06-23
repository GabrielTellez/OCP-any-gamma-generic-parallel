#include "arprec/mp_real.h"
#include "mpreal.h"
#include <gmpxx.h>

template long double cpartition::product<long double,unifunction<long double> >(int, int, unifunction<long double>*);
template mp_real cpartition::product<mp_real,unifunction<mp_real> >(int, int, unifunction<mp_real>*);
template mpfr::mpreal cpartition::product<mpfr::mpreal,unifunction<mpfr::mpreal> >(int, int, unifunction<mpfr::mpreal>*);
template mpf_class cpartition::product<mpf_class,unifunction<mpf_class> >(int, int, unifunction<mpf_class>*);

template long double cpartition::sum<long double, function<long double> >(int, int, function<long double>*, size_t);
template mp_real cpartition::sum<mp_real, function<mp_real> >(int, int, function<mp_real>*, size_t);



