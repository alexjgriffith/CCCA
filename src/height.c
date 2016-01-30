/* This file is part of CCCA,
   http://github.com/alexjgriffith/CCCA/, 
   and is Copyright (C) University of Ottawa, 2015. It is Licensed under 
   the three-clause BSD License; see LICENSE.txt.
   Author : Alexander Griffith
   Contact: griffitaj@gmail.com */

#include "height.h"


int getChromosome(char value[10])
{
  int i;
  for(i=0;i<25;i++)
    if(!strcmp(value,hg19Chrom[i].value))
      return i;
  return -1;
}

int getChromosomeShort(char value[2])
{
  int i;
  for(i=0;i<25;i++)
    if(value[0]==hg19Chrom[i].shortValue[0] && value[1]==hg19Chrom[i].shortValue[1])
      return i;
  return -1;
}

void getChromosomeValue(int order,char * value)
{
  int i;
  for(i=0;i<25;i++)
    if(order==hg19Chrom[i].order)
      {
	strcpy(value,hg19Chrom[i].value);
	return ;
      }
  strcpy(value,hg19Chrom[25].value);
}

void rankChromosomes(char ** chroms,int *length, int *out)
{
  int i;
  char string[256];
  for(i=0;i<*length;i++)
    {
      sscanf(chroms[i],"%s",string);
      out[i]=getChromosome(string);
    }
}

void valueChromosomes(char ** chroms,int *length, int *out)
{
  int i;
  for(i=0;i<*length;i++)
    {
      getChromosomeValue(out[i],chroms[i]);
    }
}

int compare(int start1, int end1,int start2, int end2)
{
  int width = (end1-start1)*2;
  int difference = (start1+end1-start2-end2);
  if(width-difference<0)
    return 0 ; //move reads forward
  else if(width+difference<0)
    return 1 ; //move peaks forward
  else
    return 2 ; //add to score and move reads forward
}

void pileup(char ** filename,char ** chro,int *start,
	    int *end,int *length,int *scores)
{
#ifdef _VERBOSE
  Rprintf("N=2\tfilename=%s\n",*filename);
#endif
#ifdef _test
  int i;
  for(i=0;i<*length;i++)
    Rprintf("%s\n",chro[i]);
#else
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  char string[1024];
  int inStart,inEnd,inChro;
  int i=0;

  int coll;
  int chrC;
  int ed1C;
  int ed2C;

  while(fgets(buffer,1024, f))
    {      
      sscanf(buffer,"%s\t%d\t%d", string,&inStart,&inEnd);
      chrC=strcmp(chro[i],string); // < 0 next peak
                                   // == 0 next compare
                                   // > 0 next reed

      while(chrC<=0){
	if(i>=*length)
	  break;
	ed1C=inStart-end[i]; //  =< 0 -- inStart < end
	                           //  > 0 next read
	ed2C=start[i]-inEnd; // =< 0-- inEnd > start
                                    // > 0 -- next peak
	if(chrC==0){
	  if(ed1C<=0 && ed2C<=0){
	    scores[i]++;
	    //coll++
	    break;
	    //i++;
	  }
	  else if (ed1C>0){
	    i++;
	  }
	  else{
	    break;
	  }
	    
	}
	else{
	  i++;
	}
	chrC=strcmp(chro[i],string);
	}
      if(i>=*length)
	break;	
    }
  fclose(f);
#endif
}

void file_length(char ** filename,int * i)
{
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  while(fgets(buffer,1024, f))
    {      
      *i=*i+1;
    }
}

void getChroms(char ** filename,char ** chroms)
{

  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  char *s;
  char chrom[256];  
  int i=0,start,end;
  int k;  
  int logic=0;
  fgets(buffer,1024, f);
  s = strtok (buffer,"\n");
  sscanf(s,"%s\t%d\t%d",chrom,&start,&end);
  strcpy(chroms[i],chrom);
  i=i+1;
  while(fgets(buffer,1024, f))
    {      
      s = strtok (buffer,"\n");
      sscanf(s,"%s\t%d\t%d",chrom,&start,&end);
      logic=0;
      for(k=0;k<i;k++)
	if( strcmp(chroms[k],chrom)==0)
	  {
	    logic=1;
	    break;
	  }
      if(logic==0)
	{
	  strcpy(chroms[i],chrom);
	  i=i+1;
	}
    
    }
}

void read_bed(char ** filename,char ** chrom,int *start, int *end)
{
  char  buffer[1024];
  FILE  * f = fopen(*filename,"r");
  char  *s;
  char string[256];
  int t;
  int i=0;
  while(fgets(buffer,1024, f))
    {      
      s = strtok (buffer,"\n");
      //sscanf(s,"%s\t%d\t%d",string,&start[i],&end[i]);
      //t=getChromosome(string);
      //getChromosomeValue(t,chrom[i]);
      sscanf(s,"%s\t%d\t%d",chrom[i],&start[i],&end[i]);
      i=i+1;
    }
  fclose(f);
}

