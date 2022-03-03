#include<algorithm>
#include<vector>
#include <cstdlib>
//contains all the type definition
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

//contains all the allocations
#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif

struct dats{
  int chr;
  double pos;
  int id;
};

bool comp(const dats& first,const dats& second){
  if(first.chr==second.chr)
    return first.pos<second.pos;
  else
    return first.chr<second.chr;
}

bool posComp(const dats& first,const dats& second){
  return first.pos<second.pos;
}

void mysort(myPars *pars,int print_info){
  try{
    if(print_info)
      flush_print("Will sort input data...                        ");
    
    std::vector<dats> data;
    // we know the endsize of the vector so lets make it this size to avoid nasty resizes
    data.reserve(pars->data->y); 
    for(int i=0;i<pars->data->y;i++){
      dats tmp = {pars->chr->array[i],pars->position->array[i],i};
      data.push_back(tmp);
    }
    //sort data
    std::sort(data.begin(),data.end(),comp);
    if(print_info>2){
      fprintf(stderr,"\n\t-> Will check uniqueness of positions:");
      
    }  
    //check for uniqueness of positions per chromosome
    for(unsigned int i=0;i<data.size()-1;i++)
      if(data[i].pos==data[i+1].pos){
	fprintf(stderr,"\t-> Non unique position: chromosome:%d position:%fMb\n ",data[i].chr,data[i+1].pos);
      }
    
    
    dArray *newPos = allocDoubleArray(pars->position->x);
    iArray *newChr = allocIntArray(pars->chr->x);
    usiMatrix *newData = allocUSIntMatrix(pars->data->x,pars->data->y);
    for(unsigned int i=0;i<data.size();i++) { 
      dats v= data[i];
      newPos->array[i] = v.pos;
      newChr->array[i] = v.chr;
      for(int s =0;s<newData->x;s++)
	newData->matrix[s][i] = pars->data->matrix[s][v.id];
      //fprintf(stderr,"chr:%d \tpos:%f \t id:%d \n",v.chr,v.pos,v.id);
    }
    killArray(pars->chr);
    killMatrix(pars->data);
    killArray(pars->position);
    
    pars->data=newData;
    pars->position=newPos;
    pars->chr=newChr;
    if(print_info)
      flush_print("\r\t-> Done sorting...                                                                 \n");
  } catch(std::exception & e){
    fprintf(stderr,"\n\t-> Allocation failed, likely cause: need more memory, will exit.\n");
    exit(0);
  } 
}
