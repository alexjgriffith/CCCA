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

void unityOutput(int * intChr, int * intSummit, int * name, int * l1, int * l2, int * peaklength, int* peakwidth, int* retChr, int * retSummit, int * retMatrix)
{
   int broken[4]={8009,8047,8474,8480};
  int i,j,k;
  int nextSummit,temp;
  for(i=0;i<(*peaklength);i++)
    {
      
      retChr[i]=intChr[l1[i]-1];
      nextSummit=0;
      j=0;
      for(k=(l1[i]-1);k<l2[i];k++)
	{
	  retMatrix[i* (*peakwidth)+name[k]-1]=1;
	  temp=nextSummit+intSummit[k];
	  nextSummit=temp;
	  j++;
	}
      if(j==0){
	Rprintf("Divide by Zero in unityOutput");
	exit(0);
      }
      retSummit[i]=nextSummit/j;
      for(k=0;k<4;k++)
	if (i== broken[k]){
	  Rprintf("chr=%d\tsummit=%d\n",retChr[i],retSummit[i]);

	}
      
    }  
}


  
// region takes an input list of form 
// value logical category categories outMatrix
// logical indicates that the value is a
// start or end value
// it will return a matrix of integers width = 2+#cat
