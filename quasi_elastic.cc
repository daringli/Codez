#include <vector>
#include <cstdio>
#include "nucleus.hh"
#include "globals.hh"
#include "lorentz_boost.hh"
#include "vector34.hh"

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp> //uniform randomness!

//Set outgoing momenta from a quasi-elastic scattering event given
//a prefragment E* and J.

std::vector<Nucleus> quasi_elastic(Nucleus projectile, Nucleus target, Nucleus cluster, double prefragment_Eex, double prefragment_J)
{
  Nucleus prefragment;
 //fragment, the participant part.
  prefragment.set_Z(projectile.Z()-cluster.Z());
  prefragment.set_N(projectile.N()-cluster.N());

  //masses
  double projectile_mass=mass_model.mass(projectile);
  double target_mass=mass_model.mass(target);
  double cluster_mass=mass_model.mass(cluster); 
  double fragment_mass=mass_model.mass(prefragment); //mass model does not include E*
  double prefragment_mass=fragment_mass+prefragment_Eex;  //prefragment mass includes E*

  prefragment.set_E(prefragment_Eex);
  prefragment.set_J(prefragment_J);
  
  //some reaction properties. Should be prefragment rather than fragment??
  double Q=-(prefragment_mass + target_mass - projectile_mass);
  double mu=prefragment_mass * target_mass/(prefragment_mass + target_mass);

  //----------END INITIALIZATION----------------
  //----------START CALCULATIONS---------------------

  //we get the initial cluster 3-momentum from "Goldhaber model".
  //essentially Gaussian distribution 
  //width sigma is not really calculated from Goldhaber here...
  double sigma=sqrt(-2*Q*mu); //as Leonid does it... not really Goldhaber!
  boost::random::normal_distribution<double> gauss(0.0, sigma);

  //initial cluster momentum in projectile frame
  Vector4 cluster_i_momentum_projectile;
 redop:
  cluster_i_momentum_projectile.v3.x=gauss(gen);
  cluster_i_momentum_projectile.v3.y=gauss(gen);
  cluster_i_momentum_projectile.v3.z=gauss(gen);
  //need off-shell mass to calculate 4-momentum. Given by 4-mom conservation
  double cluster_off_shell_mass = sqrt(pow(projectile_mass,2)+pow(prefragment_mass,2) - 2*projectile_mass *sqrt(pow(prefragment_mass,2)+cluster_i_momentum_projectile.v3*cluster_i_momentum_projectile.v3));
  //cluster energy in projectile frame as 0th component
  cluster_i_momentum_projectile.x0=sqrt(pow(cluster_off_shell_mass,2)+cluster_i_momentum_projectile.v3*cluster_i_momentum_projectile.v3);

  //printf("on-shell mass: %f\n",cluster_mass);
  //printf("off-shell mass: %f\n",cluster_off_shell_mass);


  //prefragment 3-momentum is minus cluster<=>projectile frame is ZM frame
  Vector4 prefragment_P;
  prefragment_P.v3=-cluster_i_momentum_projectile.v3;
  prefragment_P.x0=sqrt(pow(prefragment_mass,2)+prefragment_P.v3*prefragment_P.v3);

  //we now get the cluster & prefrag 4-momentum in the lab (target) frame.
  //First calculate gamma of the projectile
  double gamma_projectile=projectile.P().x0/projectile_mass; // gamma=E/m
  
  Vector4 cluster_i_momentum_lab = lorentz_boost_z(cluster_i_momentum_projectile,gamma_projectile); //z is in projectile (beam) direction.
  //We can calculate the final prefragment lab frame momentum. Set it!
  prefragment.set_P(lorentz_boost_z(prefragment_P,gamma_projectile));

  //initial target momentum in lab frame is trivial
  Vector4 target_i_momentum_lab;
  target_i_momentum_lab.v3.x=0;
  target_i_momentum_lab.v3.y=0;
  target_i_momentum_lab.v3.z=0;
  target_i_momentum_lab.x0=target_mass;
  
  
  //generate mandelstam variables: t from dsigma/dt = "uniform", and s from known quantities (invariant mass). 
  
  double s=pow(cluster_off_shell_mass,2)+pow(target_mass,2)+2*target_mass*cluster_i_momentum_lab.x0;

  //center of mass energy and p=|3-momentum| expressed in terms of invariants (s)
  double cluster_i_energy_cm = (s+pow(cluster_off_shell_mass,2)-pow(target_mass,2))/(2*sqrt(s));
  double cluster_f_energy_cm = (s+pow(cluster_mass,2)-pow(target_mass,2))/(2*sqrt(s));
  double cluster_f_p_cm = sqrt(pow(cluster_f_energy_cm,2)-pow(cluster_mass,2));
  double cluster_i_p_cm = sqrt(pow(cluster_i_energy_cm,2)-pow(cluster_off_shell_mass,2));

 if(isnan(cluster_f_p_cm))
    goto redop;
 if(isnan(cluster_i_p_cm))
    goto redop;
 if(isnan(cluster_f_energy_cm))
    goto redop;
 if(isnan(cluster_i_energy_cm))
    goto redop;

  //also need target
  double target_energy_cm = (s+pow(target_mass,2)-pow(cluster_mass,2))/(2*sqrt(s));



  double t;
 redot:
  boost::random::uniform_real_distribution<double> dist_t(0,1);  
  t=-s*dist_t(gen);
  //printf("t: %f\n",t); 
  double theta_cm = acos((t-pow(cluster_off_shell_mass,2)-pow(cluster_mass,2) + 2*cluster_i_energy_cm*cluster_f_energy_cm)/(2*cluster_f_p_cm*cluster_i_p_cm));
  //printf("theta_cm: %f\n",theta_cm);
  if(isnan(theta_cm))
    goto redot; //THIS HAPPENS WORRYINGLY OFTEN

  //randomize momentum phi direction in CM
  boost::random::uniform_real_distribution<double> dist_phi(0,2*M_PI);
  double phi=dist_phi(gen);
  Vector4 cluster_f_momentum_cm;
  cluster_f_momentum_cm.v3.x=cluster_f_p_cm*sin(theta_cm)*cos(phi);
  cluster_f_momentum_cm.v3.y=cluster_f_p_cm*sin(theta_cm)*sin(phi);
  cluster_f_momentum_cm.v3.z=cluster_f_p_cm*cos(theta_cm);
  //energy included to interface with Lorentz boost function
  cluster_f_momentum_cm.x0=cluster_f_energy_cm;
  
  Vector4 target_f_momentum_cm; 
  target_f_momentum_cm.v3= -cluster_f_momentum_cm.v3;
  //energy included to interface with Lorentz boost function
  target_f_momentum_cm.x0= target_energy_cm;


  //boost back to lab frame
  Vector3 beta_cm_to_lab=-cluster_i_momentum_lab.v3*(1/(cluster_i_momentum_lab.x0 + target_mass)); 
  Vector4 cluster_f_momentum=lorentz_boost(cluster_f_momentum_cm,beta_cm_to_lab);
  Vector4 target_f_momentum=lorentz_boost(target_f_momentum_cm,beta_cm_to_lab);

  cluster.set_P(cluster_f_momentum);
  target.set_P(target_f_momentum);
  cluster.set_E(0.0);
  target.set_E(0.0);
  std::vector<Nucleus> ret(3);
  ret[0]=prefragment;
  ret[1]=cluster;
  ret[2]=target;

  return ret;

}
