#include "region.h"

int addToBuffer(int logical, int category, int *buffer)
{
  if(logical==0) 
    buffer[category]=1;
  else 
    buffer[category]=0;
  return 0;
}

void region(int * value, int * logical , int * category,int * buffer, int * outMatrix,int * length,int * width)
{
  int start,end,i,j;
  start=value[0];
  addToBuffer(logical[0],category[0],buffer);
  for(i=1;i<*length;i++)
    {
      end=value[i];
      outMatrix[i*(*width+2)+0]=start;
      outMatrix[(i*(*width+2)+1)]=end;
      for(j=0;j<*width;j++)
	outMatrix[(i*(*width+2)+j+2)]=buffer[j];      
      addToBuffer(logical[i],category[i],buffer);
      start=end;
    }
}

  
// region takes an input list of form 
// value logical category categories outMatrix
// logical indicates that the value is a
// start or end value
// it will return a matrix of integers width = 2+#cat
