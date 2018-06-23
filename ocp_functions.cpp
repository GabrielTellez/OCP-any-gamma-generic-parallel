#include "ocp_functions.h"

template <>
long double my_pi<long double>()
{
  return M_PI;
}
template <>
mp_real my_pi<mp_real>()
{
  return mp_real::_pi;
}
