#ifndef POTENTIAL__HH
#define POTENTIAL__HH
class Potential
{
public:
  virtual double value(double r) const=0;
};
#endif
