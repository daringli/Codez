#include "file_names.hh"
#include <cstring>

char* outfile_name(int Z, int N, double E, double J, char* ending)
{
  static char filename[40]="ZxxxNxxxExxxjxxx";
  filename[1]=(char)((int)'0'+Z/100);
  filename[2]=(char)((int)'0'+(Z-100*(Z/100))/10);
  filename[3]=(char)((int)'0'+Z%10);
  filename[5]=(char)((int)'0'+N/100);
  filename[6]=(char)((int)'0'+(N-100*(N/100))/10);
  filename[7]=(char)((int)'0'+N%10);
  int iE=int(E);
  filename[9]=(char)((int)'0'+iE/100);
  filename[10]=(char)((int)'0'+(iE-100*(iE/100))/10);
  filename[11]=(char)((int)'0'+iE%10);

  filename[13]=(char)((int)'0'+((int)J)/10);
  filename[14]=(char)((int)'0'+((int)J)%10);
  filename[15]=(char)((int)'0'+((int)(10*J))%10);
  filename[16]='\0';

  //problem, can't just strcat since filename may be full of 
  //other things from previous iterations.
  strcat(filename,ending);
  return &filename[0];
}
