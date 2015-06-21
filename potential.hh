#ifndef POTENTIAL__HH
#define POTENTIAL__HH
class Potential
{
public:
  virtual double value(double r) const=0;
  virtual double first_derivative(double r) const=0;
  virtual double second_derivative(double r) const=0;

};
#endif
