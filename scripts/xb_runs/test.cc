#include <vector>
#include <cstdio>
int main(){
  int x;
  int size=3;
  for(int i=0;i<size;i++){
    *(&x+i)=i; //x={1,2,3,...}
    printf("x[%i]=%i\n",i,*(&x+i));
  }
  std::vector<int> vx(&x,&x+3);
}
