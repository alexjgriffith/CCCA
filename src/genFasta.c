#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <R_ext/Error.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
//#define _V
#define _M
#define _SAVE


void generateFasta(char **sourceFasta,char **indexFasta,char ** chrs, int * starts, int * ends,int *length, int * width,char ** fastadata);

int generateIndex(char ** sourceFasta,char ** indexFasta);

int generateIndex(char ** sourceFasta,char ** indexFasta)
{
  FILE * file;
  FILE * write;
  int i;
  char buffer[10];
  char chro[10];
  int k;
  int n=0;
  file=fopen(*sourceFasta,"r");
  write=fopen(*indexFasta,"w");
  if(file==NULL) 
    error("File, %s, does not exist",*sourceFasta);
  if(write == NULL)
    error("File, %s, does not exist",*indexFasta);
  while (fread(buffer,1,1,file))
    {
      if(n%100000==0)
	R_CheckUserInterrupt();
      if (buffer[0]=='>')
	{
	  for(k=0;k<10;k++)
	    chro[k]='\0';
	  k=0;  
	    while(fread(buffer,1,1,file))
	      {
		if(buffer[0]=='\n' || k==9)
		  break;
		chro[k]=buffer[0];
		k++;
	      }
	    fprintf(write,"%s\t%lu\n",chro,ftell(file));
	}
      n++;
    }
  fclose(write);
  fclose(file);
  return 0;
}



void generateFasta(char **sourceFasta,char **indexFasta,char ** chrs, int * starts, int * ends,int *length, int * width,char ** fastadata)
{
#ifdef _V
  Rprintf("_V> Verbose Flag (_V) in genFasta.c is defined\n");
#endif

  FILE * file;
  FILE * fasta;
  char t2[2048];
  char *buffer =(char * )R_alloc(256,sizeof(char));
  char tempstr [2048];
  int k;
  int i;
  int j;
  int chroLength;
  char caps[100][100];
  unsigned long temp[100];
  unsigned long rats,line;
  char chro[100];
  unsigned long start;

#ifdef _M
  Rprintf("_M> Message Flag (_M) in genFasta.c is defined\n");
  int percents[100];
  int pcount=0;
  int pmax=10;
  if(*length<100)
    pmax=1;
  else if(*length>100000)
    pmax=100;
  percents[0]=floor(*length)/pmax;
  for(i=1;i<pmax;i++)
    percents[i]=percents[i-1]+floor(*length)/pmax  ;
#endif

#ifdef _V
  Rprintf("_V> Break-1\n",temp[k]);
#endif
  file=fopen(*indexFasta,"r"); 
  if(file==NULL){
#ifdef _V
    Rprintf("file=NULL\n");
    Rprintf("Could not find the index file: %s\ngenerating %s now.\n",indexFasta,indexFasta);
#endif
    warning("Could not open index file: %s\n",*indexFasta);
    generateIndex(sourceFasta,indexFasta);
    file=fopen(*indexFasta,"r"); 
  }
  fasta=fopen(*sourceFasta,"r");     
#ifdef _V
  Rprintf("_V> Break-2\tindexFile=%s\n",*indexFasta);
#endif

  k=0;
  while (fgets(buffer,100,file))
    {
      sscanf(buffer, "%s\t%lu",caps[k],&temp[k]);
#ifdef _V
      Rprintf("_V> %s\t%lu\n",caps,temp[k]);
#endif
      k++;
    }
#ifdef _V
  Rprintf("_V> Break-3\n",temp[k]);
#endif

  fclose(file);
  chroLength=k;
#ifdef _V
  Rprintf("_V> Break-4\n",temp[k]);
#endif
#ifdef _SAVE
  FILE * outfile=fopen(*fastadata,"w");
#endif
  for(j=0;j<(*length);j++)
    {
#ifdef _M
      if(j==percents[pcount])
	{
      Rprintf("_M> Progress = %d\%\n",(100*(1+pcount))/pmax);
      if(pcount<pmax)
	pcount++;
	}
#endif
      if(j%1000==0)
	R_CheckUserInterrupt();
      start=(unsigned long) (starts[j]+ends[j]-(*width))/2;
      for (i=0;i<chroLength;i++)
	{
	  sscanf(chrs[j],"%s",chro);
	  if (!strcmp(caps[i],chro))
	    {
	      line=(unsigned long)start/50;	      
	      rats=start+temp[i]+line;
	      fseek(fasta,rats,SEEK_SET);
	      for(k=0;k<(*width);k++)
		{
		  tempstr[k]='\0';
		}
	      for(k=0;k<(*width);k++)
		{
		  if(fread(buffer,1,1,fasta)){
		  
		    if (buffer[0]=='\n')
		      k--;
		    else     
		      tempstr[k]=buffer[0];
		  }		  
		  else		   
		    break;
		}

	      tempstr[k+1]='\0';
	      sscanf(tempstr,"%s",t2);
#ifndef _SAVE
	      strcpy(fastadata[j],t2);
#endif
#ifdef _SAVE
	      fprintf(outfile,">%s:%d-%d\n%s\n",chrs[j],starts[j],ends[j],t2);
#endif
	      break;
	    }
	}
    }
  fclose(fasta);
#ifdef _SAVE
  fclose(outfile);
#endif
#ifdef _V
  Rprintf("_V> ending\n");
#endif

}



