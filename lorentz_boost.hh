#ifndef LORENTZ_BOOST_Z__HH
#define LORENTZ_BOOST_Z__HH

#include "vector34.hh"

Vector4 lorentz_boost(Vector4 v, Vector3 beta);

Vector4 lorentz_boost_z(Vector4 v, double gamma);
#endif
