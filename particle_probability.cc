#include <cstdio>

#include "particle_probability.hh" 
#include "decay.hh"
#include "decay_name.hh"

std::vector<Particle_probability> particle_width(std::vector<std::pair<double,Decay> > &table)
{
  //get total decay width to normalize with
  int N=table.size();
  double Gamma_tot=table[N-1].first;
  
  std::vector<Particle_probability> pps(0);

  //get partial width of different particle channels in the order the appear in the decay table
  struct {int N; int Z;   bool operator!=(const Decay& a) const{return (N != a.N) && (Z!=a.Z);}} channel;

  //initial channel, since table[0] is padding.
  channel.N=table[1].second.N;
  channel.Z=table[1].second.Z;  
  double previous_cum_width=0; //cummulative decay width to previous decay.
  Particle_probability this_pp;
  for(int i=1;i<N;i++){
    if(channel != table[i].second){
      //next channel starts here!
      this_pp.prob=(table[i-1].first-previous_cum_width)/Gamma_tot;
      previous_cum_width=table[i-1].first;    
      this_pp.name=decay_name(table[i-1].second);
      pps.push_back(this_pp);
      channel.N=table[i].second.N;
      channel.Z=table[i].second.Z;  
    }
  }
  //the final channel
  this_pp.prob=(table[N-1].first-previous_cum_width)/Gamma_tot;
  this_pp.name=decay_name(table[N-1].second);
  pps.push_back(this_pp);

  return pps;
}

void print_table(std::vector<std::pair<double,Decay> > &table, Nucleus n)
{
  //get total decay width to normalize with
  int N=table.size();
  puts("#Zev\tNev\tJf\tEi-Ef\tl\tp");//for some reason, the program segfaults without this
  double Gamma_tot=table[N-1].first;
  double Ei=n.E();
  double previous_cum_width=0; //cummulative decay width to previous decay.
  for(int i=1;i<N;i++){
    if(i>1){
      previous_cum_width=table[i-1].first;    
    }
    printf("%i\t%i\t%f\t%f\t%i\t%e\n",table[i].second.Z,table[i].second.N,table[i].second.j*0.5,Ei-table[i].second.E,table[i].second.l,table[i].first-previous_cum_width);
  }
}
