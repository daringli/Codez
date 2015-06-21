#include "integrate.hh"
#include <cmath>

double integrate(Potential const & pot, double E, double ri, double ro)
{
  //following MyCodex, we split up the interval in parts and use
  //Guassian quadrature with n=7 on each subinterval.
  //Integrate sqrt(V-E)
  int intervals=2;
  int n=7;
  int index; //since a[-j]=a[j], we use index=abs(j).
  int values=n/2+n%2; //unique tabulated values

  double subintervalLength=(ro-ri)/intervals; 
  double r1=ri-subintervalLength; 
  double r2;
  double center;

  double a[4]={0.417959183673469,  0.381830050505119, 0.279705391489277, 0.129484966168870};
  double x[4]={0.0, 0.405845151377397,0.741531185599394,0.949107912342759};

  double integral=0;
  double r;
  for(int i=1;i<=intervals;i++){
    r1=r1+subintervalLength;
    r2=r1+subintervalLength;
    center=r1+0.5*subintervalLength;
      
    for(int j=-values; j<=values;j++){
      index=abs(values);
      r=x[index]*0.5*subintervalLength+center;
      integral+=0.5*subintervalLength*(sqrt(pot.value(r)-E));
    }
  }
  return integral;
}
