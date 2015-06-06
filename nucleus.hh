#ifndef NUCLEUS__HH
#define NUCLEUS__HH

#include <cmath>       //pow, exp, M_PI, sqrt
#include "vector34.hh"

class Nucleus
{
private:
  int m_Z;
  int m_N;
  double m_E;
  double m_J;
  double m_gs_beta; //beta for ground state. Gives quadrupole deformation. R=R_0(1+beta2*Y^2_0(theta))
  double m_gs_deformation_energy;
  Vector4 m_P;
public:

  Nucleus(int z, int n){
   //constructor that sets values that tells us that the user didn't specify a beta and fissility parameter. -1 perhaps unsuitable beta sentry value. (beta can take on negative values, but Möller only gets beta~-0.1 at worst.)
    m_Z=z;
    m_N=n;    
    m_gs_beta=-1;  
  }
 Nucleus(){
   
    m_gs_beta=-1;    
  }
  Nucleus(int z, int n,double bGS){
    //constructor to let user specify deformations
    m_Z=z;
    m_N=n;   
    m_gs_beta=bGS;    
  }
  Nucleus(double bGS){
    //constructor to let user specify deformations
    m_gs_beta=bGS;    
  }

  //getters
  int A() const{return m_N+m_Z;}
  int Z() const{return m_Z;}
  int N() const{return m_N;}
  double E()const{return m_E;}
  double J() const{return m_J;}
  double gs_beta() const{return m_gs_beta;}
  double gs_alpha() const{int l=2; return sqrt((2.0*(double)l+1.0)/M_PI)*m_gs_beta;}
 //From Spherical Harmonics to Legendre Polynomial normalization.
  double gs_deformation_energy() const{return m_gs_deformation_energy;}
  Vector4 P() const{return m_P;}

  //setters
  void set_Z(int Zs){m_Z=Zs;}
  void set_N(int Ns){m_N=Ns;}
  void set_E(double E){m_E=E;}
  void set_J(double Js){m_J=Js;} 
  void set_J(int Js){set_J((double)Js);} 
  void set_gs_deformation_energy(double Edef){m_gs_deformation_energy=Edef;}
  void set_P(Vector4 v){m_P = v;}
};

#endif
