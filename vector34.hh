#ifndef vector34__hh
#define vector34__hh

struct Vector3
{
  double x;
  double y;
  double z;

  Vector3 operator=(const Vector3& a) {
        x=a.x;
        y=a.y;
	z=a.z;
        return a;
  }

  Vector3 operator+(const Vector3& a){
    Vector3 v;
    v.x = this->x + a.x;
    v.y = this->y + a.y;
    v.z = this->z + a.z;
    return v;
  }
  Vector3 operator-(const Vector3& a){
    Vector3 v;
    v.x = this->x - a.x;
    v.y = this->y - a.y;
    v.z = this->z - a.z;
    return v;
  }

  //scalar product
  double operator*(const Vector3& a){
    return this->x*a.x + this->y*a.y + this->z*a.z;
  }

  //scale with integer
  Vector3 operator*(const double a){
    Vector3 v;
    v.x = this->x*a;
    v.y = this->y*a;
    v.z = this->z*a;
    return v;
  }
  Vector3 operator-() const {
      Vector3 v;
      v.x = -x;
      v.y = -y;
      v.z = -z;
      return v;
   }
};


struct Vector4
{
  double x0;
  Vector3 v3;

  Vector4 operator=(const Vector4& a) {
        x0=a.x0;
        v3=a.v3;
        return a;
  }

  Vector4 operator+(const Vector4& a){
    Vector4 v;
    v.x0 = this->x0 + a.x0;
    v.v3 = this->v3 + a.v3;
    return v;
  }
  Vector4 operator-(const Vector4& a){
    Vector4 v;
    v.x0 = this->x0 - a.x0;
    v.v3 = this->v3 - a.v3;
    return v;
  }
 

  //scalar product
  double operator*(const Vector4& a){
    return this->x0*a.x0-this->v3*a.v3; 
  }
    //scale with integer
  Vector4 operator*(const double a){
    Vector4 v;
    v.x0 = this->x0*a;
    v.v3 = this->v3*a;
    return v;
  }
};

#endif
