#include "vector34.hh"
#include <cmath>

Vector4 lorentz_boost(Vector4 v, Vector3 beta)
{
  double b=sqrt(beta*beta);
  double gamma=1/sqrt(1-b*b);
  Vector4 ret;

  ret.x0=gamma*(v.x0-beta*v.v3);
  ret.v3.x=-gamma*beta.x*v.x0 +(1+(gamma-1)*beta.x*beta.x*pow(1/b,2))*v.v3.x +((gamma-1)*beta.x*beta.y*pow(1/b,2))*v.v3.y +((gamma-1)*beta.x*beta.z*pow(1/b,2))*v.v3.z;
  ret.v3.y=-gamma*beta.y*v.x0 +(1+(gamma-1)*beta.y*beta.y*pow(1/b,2))*v.v3.y +((gamma-1)*beta.y*beta.x*pow(1/b,2))*v.v3.x +((gamma-1)*beta.y*beta.z*pow(1/b,2))*v.v3.z;
  ret.v3.z=-gamma*beta.z*v.x0 +(1+(gamma-1)*beta.z*beta.z*pow(1/b,2))*v.v3.z +((gamma-1)*beta.z*beta.x*pow(1/b,2))*v.v3.x +((gamma-1)*beta.z*beta.y*pow(1/b,2))*v.v3.y;
  return ret;
}

Vector4 lorentz_boost_z(Vector4 v, double gamma)
{
  double beta=-sqrt(1-1/pow(gamma,2));
  Vector4 ret;
  ret.v3=v.v3;
  ret.x0=gamma*(v.x0-beta*v.v3.z);
  ret.v3.z=gamma*(v.v3.z-beta*v.x0);
  return ret;
}
