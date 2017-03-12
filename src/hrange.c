#include "hrange.h"

void peakDensity(char ** filename,char ** chro,int *start,
	    int *end,int *length,int *scoresOut){
  int j;
  int i;
  peak *temp;
  int * lengths= malloc(*length*sizeof(int));
  peak ** scores= malloc(*length*sizeof(int*));
  int ** collect=malloc(*length*sizeof(int*));

  buildPeaks(chro,start,end,*length,&temp);
  //Rprintf("out\n");
  //Rprintf("%d\n",start[0]);
  //Rprintf("%d\n",temp[0].start);
  for(i=0;i<*length;i++){
    lengths[i]=0;
    scores[i]=malloc(100*sizeof(peak));    
    collect[i]=malloc((temp[i].end-temp[i].start) *sizeof(peak));
  }
  //printf("Allocated\n");

  rheight(*filename,temp,*length,&scores,&lengths);

  convertHeights(temp,*length, scores, lengths, &collect);
  //printConvertedHeights( collect, temp, *length);
  int width=temp[0].end-temp[0].start;
  for(i=0;i<*length;i++){    
    for(j=0;j<width;j++)
	scoresOut[i*width+j]=collect[i][j];
  }
}

void buildPeaks(char **chr, int * start, int * end , int length, peak ** temp)
{
  int i,k;
  (*temp)=malloc(sizeof(peak)*length);

  for(i=0;i<length;i++){
    k=0;
    while(chr[i][k]!='\0'){
      (*temp)[i].chr[k]=chr[i][k];
      k++;
    }
    (*temp)[i].chr[k]='\0';

    (*temp)[i].start=start[i];

    (*temp)[i].end=end[i];
  }
  //Rprintf("%s:%d-%d\n",(*temp)[0].chr,(*temp)[0].start,(*temp)[0].end);
}

void rheight(char * filename,peak * peaks,int length, peak *** scores,int ** heights){
  char  buffer[1024];
  FILE  * f = fopen(filename,"r");
  char string[1024];
  int inStart,inEnd;
  int i=0;
  int j=0;
  int count=0;
  int coll;
  int chrC;
  int ed1C;
  int ed2C;  
  peak holding[100];
  peak temp;
  int expanded=0;
  while(fgets(buffer,1024, f))
    {   
      //printf("%d\t$d\n",i,count);      
      sscanf(buffer,"%s\t%d\t%d", string,&inStart,&inEnd);
      chrC=strcmp(peaks[i].chr,string);
      // if there are more than 100 peaks assosiated expand memory
      if(count>98 && expanded==0){
	//printf("expanding memory\n");
	for(j=0;j<100;j++)
	  holding[j]=(*scores)[i][j];
	free((*scores)[i]);
	(*scores)[i]=malloc(2000*sizeof(peak));
	for(j=0;j<100;j++)
	  (*scores)[i][j]=holding[j];
	//printf("memory expanded\n");
	expanded=1;
	}	

      while(chrC<=0){

	if(i>=length)
	  break;
	ed1C=inStart-peaks[i].end; //  =< 0 -- inStart < end
	                           //  > 0 next read
	ed2C=peaks[i].start-inEnd; // =< 0-- inEnd > start
                                    // > 0 -- next peak
	if(chrC==0){	  
	  if(ed1C<=0 && ed2C<=0){
	    //printf("adding read to peak\n");
	    strcpy(temp.chr,peaks[i].chr);
	    temp.start=inStart;
	    temp.end=inEnd;
	    (*scores)[i][count]=temp;
	    (*heights)[i]++;
	    //printf("peak added to read\n");
	     count++;
	    break;
	  }
	  else if (ed1C>0){
	    expanded=0;
	    count=0;
	    i++;
	  }
	  else{
	    break;
	  }
	    
	}
	else{
	  expanded=0;
	  count=0;
	  i++;
	}
	chrC=strcmp(peaks[i].chr,string);
	}
      if(i>=length)
	break;	
    }

  fclose(f);

}


void convertHeights(peak * temp, int  length, peak ** scores, int * lengths,int *** collectIn){
  int height=0;
  int ht;
  int j;
  int k;
  int i;
  for(i=0;i<length;i++){
    for(k=-(temp[i].end-temp[i].start);k<2*(temp[i].end-temp[i].start);k++){
      for(j=0;j<lengths[i];j++){
	if((temp[i].start+k)==scores[i][j].start){
	  ht=height;
	  height=ht+1;
	}
	if((temp[i].start+k)==scores[i][j].end){
	  ht=height;
	  height=ht-1;	
	}
      }
      if(height<0)
	error("Height cannot be less than 0\n");
      if(k>=0 && k< (temp[i].end-temp[i].start)){
	(*collectIn)[i][k]=height;

      }
    }
  }
}

void printConvertedHeights(int ** collect,peak * temp, int length){
  int i;
  int j;
  for(i=0;i<length;i++){
    Rprintf("%s:%d-%d",temp[i].chr,temp[i].start,temp[i].end);
    for(j=0;j<(temp[i].end-temp[i].start);j++){
	Rprintf("%i ",collect[i][j]);      
    }
    Rprintf("\n");
  }
}

