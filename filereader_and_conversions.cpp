#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif
#define PLINK_POS_SCALING 1000000 //used for scaling the position from a plink format
#define MAX_ELEMS_PER_LINE 5000000

 
#include <string>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "filereader_and_conversions.h"


void fileError(std::string file){
  printf("\t-> Problems opening file: %s\t right directory?\n",file.c_str());
}


void getOptions(const char *fileName, functionPars *pars){
  float options[29];		//array of option arguuments
  std::ifstream exsamplefile(fileName);	//opens option file
  //  char *tmpChar;
  if (!exsamplefile.is_open ()){
    printf ("\nError opening options file\n");
    exit(0);
  }
  int tal = 0;			// read option file
  char buffer[256];
  while (!exsamplefile.eof ()) {
    exsamplefile.getline (buffer, 1000);
    char *tokenized=buffer;
    char *buffer_ar = strtok(tokenized,"\t"); //changed from strsep in v. 0.82
    
    if(tal<27){
      //printf("tal:%d\tbuffer:%s\n",tal,buffer_ar);
      options[tal] = atof (buffer_ar);
      tal++;
      continue;
    }
    else{
      // there exist more lines in options file let's check
      if(buffer_ar!=NULL){//The lines are not empty, lets discard them but tell user
	printf("tal:%d\t op:%f \tlen=%d\n",tal,options[tal],(int)strlen(buffer_ar));
	printf("\t-> To many lines in options file, rest will be discarded\n");
	break;
      }
      //lines are empy just keep on going
    }
  }

  exsamplefile.close ();

  pars->doAllPairs = (int) options[0];
  pars->pair->array[0]=(int) options[1];
  pars->pair->array[1]=(int) options[2];
  pars->double_recom=(int) options[3];
  pars->LD = ((int) options[4]==0?1:0 );//ld=1=rsq, ld=0=D
  pars->min=(double) options[5];
  pars->alim->array[0]= options[6];
  pars->alim->array[1]= options[7];
  pars->doPar=(int) options[8];
  pars->par->array[0] = options[9];
  pars->par->array[1] = options[10];
  pars->par->array[2] = options[11];
  pars->ld_adj= (int) options[12];
  pars->epsilon=  options[13];
  pars->back= (int) options[14];
  pars->doPrune = (int) options[15];
  pars->prune_val =  options[16];
  pars->fixA= (int) options[17];
  pars->fixA_val=  options[18];
  pars->fixK= (int) options[19];
  pars->fixK_val= options[20];
  pars->calcA= (int) options[21];
  pars->phi =  options[22];
  pars->convTol = options[23];
  pars->timesToConverge=(int) options[24];
  pars->timesToRun= (int) options[25];
  pars->back2 = (int) options[26];
  
  
}



/// Checking for file existence, using stat.h.
int fexists(std::string str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str.c_str(), &buffer )==0 ); /// @return Function returns 1 if file exists.
}

///Converts a integer type to string.
std::string to_string(int a){///@param a An integer to convert.
  std::string str;
  std::stringstream ss;
  ss<<a;
  if(!(ss>>str)){
    printf("Error converting :%d will exit\n",a);
    exit(0);
  }
  ss>>str;
  return str; ///@return Returns the string value.
}

///Converts a string to integer type.
int to_int(std::string a){///@param a a String to convert.
  int str;
  float check;
  std::stringstream ss;
  std::stringstream ss2;
  if(!a.compare("NA")||!a.compare("Na")||!a.compare("na")){
    return 0;
  }
  ss<<a;
  if(!(ss>>str)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }

  //now check if numberis really an int
  ss2<<a;
  if(!(ss2>>check)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }

  if(check!=str){
    printf("String you want as in integer is not really an integer: \"%s\"\n",a.c_str());
    puts("Just a warning");
  }
  //check values
  /* this was commented in v. 0.983 since the chromosome file converter uses it.
  if(str<0||str>3){
    printf("Value in genotype is out of bounds, should be {0...3} value is:%d\n",str);
    exit(0);
  }
  */
  return str;///@return Function returns the integer value.
}

///Convert string to float.
float to_float(std::string a){///@param a is the string to convert.
  float str;
  std::stringstream ss;
  ss<<a;
  if(!(ss>>str)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }
  ss>>str;
  return str;///@return The float value.
}

///convert string to double
double to_double(std::string a){///@param a is the string to convert.
  double str;
  std::stringstream ss;
  ss<<a;
  if(!(ss>>str)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }
  ss>>str;
  return str;///@return The double value.
}


/*
  Mother-of-all string tokenizers
  1. args = string to be splittet
  2. args = a vector<string> reference that will contain the splittet string
  3. args = delimiter list in a string like 
  
  get_lexemes(STRING_TO_SPLIT,RESULTS_STRING_VECTOR,DELIM=" \t:,;")
  this will split on whitespace,tab,colons,comma and semicolons
  return value will be number of elems tokenized
*/
///Mother-of-all string tokenizers! Function will split a string from a given set of delimiters.
int get_lexemes(std::string str,std::vector<std::string> &lexems,const std::string delims){
  ///@param str The string to be splittet.
  ///@param lexems A reference to a list of strings in which the splittet string will be put as seperate elements. A queue.
  ///@param delims A string containing all the chars that the function will split by.
  char *token = strtok( const_cast<char*>( str.c_str() ),delims.c_str());
  int numToks = 0;
  while ( token != NULL ){
    numToks++;
    lexems.push_back( token);
    token = strtok( NULL,delims.c_str() );
  }
  return numToks;///@return The number of strings that has been splittet.
}


void typecast_stringarray_to_int_matrix(const std::vector<std::string> &tokens,iMatrix *mat){
  int start=0;
  for (int i=0;i<mat->x;i++)
    for(int j=0;j<mat->y;j++){
      mat->matrix[i][j] = to_int( tokens.at(start++));

    }
}

iArray *typecast_stringarray_to_iArray(const std::vector<std::string> &tokens){
  iArray *ary = allocIntArray(tokens.size());
  for (unsigned int i=0;i<tokens.size();i++)
    ary->array[i] = to_int( tokens.at(i));
  return ary;
}

dArray *typecast_stringarray_to_dArray(const std::vector<std::string> &tokens){
  dArray *ary = allocDoubleArray(tokens.size());
  for (unsigned int i=0;i<tokens.size();i++)
    ary->array[i] = to_double( tokens.at(i));
  return ary;
}


std::vector<std::string> readarray(std::string filename,const std::string delims=",;: \t"){
  ///@param filename The name of the file to open.
  ///@param delims A string of delimeters to split by.
  std::vector<std::string> tokens;

  const int SIZE = MAX_ELEMS_PER_LINE;
  char *buffer = new char[SIZE];
  std::ifstream pFile (filename.c_str(),std::ios::in);
  
  if(!pFile){
    fileError(filename);
    exit(0);
  }
  
  std::string tmp_string;
  //we are reading an array, so just keep reading until no more lines
  while(!pFile.eof()){
    pFile.getline(buffer,SIZE);
    tmp_string = std::string(buffer);
    get_lexemes(tmp_string,tokens,delims);
    
  }
  pFile.close();
  delete [] buffer;
  ///now we have a token array of string coerce the types now
  //typecast_stringarray(tokens);
  return tokens;
}



iMatrix *readmatrix_filty_memory(std::string filename,const std::string delim=",;: \t"){
  ///@param filename A filename to read.@param delim A string of delimiters.
  std::vector<std::string> tokens;

  const int SIZE = MAX_ELEMS_PER_LINE;//defined in conf.h
  char buffer[SIZE];
 
  std::ifstream pFile (filename.c_str(),std::ios::in);
 
  
  if(!pFile){
    fileError(filename);
    exit(0);
  }

  
  std::string tmp_string;
  int doFirstRow =1;
  int itemsInFirstRow=0;
  int numRows =0;
  while(!pFile.eof()){
     pFile.getline(buffer,SIZE);
    tmp_string = std::string(buffer);
    if(doFirstRow){
      //if file has a emptystart line
      itemsInFirstRow = get_lexemes(tmp_string,tokens,delim);
      
      if (itemsInFirstRow==0)
	continue;
      // printf("items in first rwo:%d\n",itemsInFirstRow);
      doFirstRow=0;
      numRows++;
    }
    else{
      int nItems = get_lexemes(tmp_string,tokens,delim);
	//if line is empty
	if(nItems==0)
	  continue;
	numRows++;
	if(nItems!=itemsInFirstRow){
	  printf("row length mismatch at line:%d numitems is:%d shouldn't be:%d\t will exit\n",numRows,itemsInFirstRow,nItems);
	  exit(0);
	}
      }
  }
  flush_print("\r\t-> File has been read in now, will now typecheck...                              ");
  iMatrix *data_ = allocIntMatrix(numRows,itemsInFirstRow);
  //now we have a token array of string coerce the types now
  typecast_stringarray_to_int_matrix(tokens,data_);
  //copy(tokens.begin(), tokens.end(), ostream_iterator<string>(cout, ", "));
  return data_; 
}




iMatrix *readmatrix(std::string filename,const std::string delim=",;: \t"){
  ///@param filename A filename to read.@param delim A string of delimiters.
  if(0){
    printf("\t-> will try to open postfile: \"%s\" ... \n",filename.c_str());
    fflush(stdout);
  }
  


  const int SIZE = MAX_ELEMS_PER_LINE;
  char buffer[SIZE];
 
  std::ifstream pFile (filename.c_str(),std::ios::in);
 
  
  if(!pFile){
    printf("\t-> will try to open postfile: \"%s\" ... \n",filename.c_str());
    //    std::cout <<"Problems opening file" <<filename<<std::endl;
    exit(0);
  }
  
  
  int doFirstRow =1;
  int itemsInFirstRow=0;
  int numRows =0;
  flush_print("Checking consistency of file...");
  
  while(!pFile.eof()){
    
    std::vector<std::string> tokens;
    std::string tmp_string;

    pFile.getline(buffer,SIZE);
    tmp_string = std::string(buffer);
    
    if(doFirstRow){
      //if file has a emptystart line
      itemsInFirstRow = get_lexemes(tmp_string,tokens,delim);
      
      if (itemsInFirstRow==0)
	continue;
      // printf("items in first rwo:%d\n",itemsInFirstRow);
      doFirstRow=0;
      numRows++;
    }
    else{
      int nItems = get_lexemes(tmp_string,tokens,delim);
	//if line is empty
	if(nItems==0)
	  continue;
	numRows++;
	if(nItems!=itemsInFirstRow){
	  printf("row length mismatch at line:%d numitems is:%d shouldn't be:%d\t will exit\n",numRows,itemsInFirstRow,nItems);
	  exit(0);
	}
      }
    if ((numRows%20 )==0){
      printf("\r\t-> Checking consistency of file: (checking number of items at line: %d )",numRows);
      fflush(stdout);
    }
  }
  pFile.close();

  fflush(stdout);

  
  iMatrix *mat = allocIntMatrix(numRows,itemsInFirstRow);

  numRows = 0;

  std::ifstream pFile2 (filename.c_str(),std::ios::in);
  while(!pFile2.eof()){
    if ((numRows%5 )==0){
      printf("\r\t-> Checking consistency of file: (Now reading in data at line: %d/%d   ) ",numRows,mat->x);
      fflush(stdout);
    }
    
    std::vector<std::string> tokens;
    std::string tmp_string;

    pFile2.getline(buffer,SIZE);
    tmp_string = std::string(buffer);
   
    int itemsInRow = get_lexemes(tmp_string,tokens,delim);
    if (itemsInRow==0)
	continue;
    for(unsigned int i=0; i <tokens.size();i++){
      mat->matrix[numRows][i] = to_int(tokens[i]);
      if(mat->matrix[numRows][i]<0||mat->matrix[numRows][i]>3){
	printf("\n\t-> Error in genotype data: (%d,%d)=%d, value should be between 0 and 3.\n",numRows,i,mat->matrix[numRows][i]);
	exit(0);
      }
    }

    numRows++;
  }
  pFile2.close();
  //copy(tokens.begin(), tokens.end(), ostream_iterator<string>(cout, ", "));
  
  fflush(stdout);
  return mat;
}





iMatrix *getData(const char * pName){
  return readmatrix(std::string(pName),"\t ;,:");
}


iArray *getChr(const char *pName){
  return typecast_stringarray_to_iArray(readarray(std::string(pName)));

}



int getNext(iArray *list,int i){
  if(i>=list->x)
    return -1;
  for(int j=i;j<list->x;j++){
    // printf("i:%d\tj:%d\tiar:%d\n",i,j,list->array[i]);
    if(list->array[j]==1)
      return j;
  }
  return -1;
}


std::vector<iArray*> getTestIndividuals(const char *pName){
  iArray *testList = typecast_stringarray_to_iArray(readarray(std::string(pName)));
  
  int numTrues=0;
  for(int i=0;i<testList->x;i++)
    if(testList->array[i])
      numTrues++;
  //printf("numTrues:%d\n",numTrues);
  //testList->print("testlist",NULL);
  std::vector<iArray*> returnVal;
  int pos1 = getNext(testList,0);
  while(pos1!=-1){
    int pos2 = getNext(testList,pos1+1);
    while(pos2!=-1){
      if(pos2!=-1){
	//printf("p1:%d\tp2:%d\n",pos1,pos2);
	iArray *pr = allocIntArray(2);
	pr->array[0]=pos1;pr->array[1]=pos2;
	returnVal.push_back(pr);
      }
      pos2 = getNext(testList,pos2+1);
    }
    pos1 = getNext(testList,pos1+1);
  }
  killArray(testList);
  return returnVal;

}



dArray *getPos(const char *pName){
  return (typecast_stringarray_to_dArray(readarray(std::string(pName))));
}


iMatrix *bed_to_iMatrix(const char* file, int nrow,int ncol) {

  //  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  //0,1,2: 3 is missing 
  const unsigned char recode[4] = { '\x02','\x01', '\x03', '\x00'};
  const unsigned char mask = '\x03';


  FILE *in = fopen(file, "r");
  if (!in){
    printf("Couln't open input file: %s", file);
    exit(0);
  }
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3){
    printf("Failed to read first 3 bytes");
    exit(0);
  }
  if (start[0]!='\x6C' || start[1]!='\x1B'){
    printf("Input file does not appear to be a .bed file (%X, %X)", 
	   start[0], start[1]);
    exit(0);
  }
  /* Create output object */
  
  
  iMatrix *returnMat = allocIntMatrix(nrow,ncol);
  unsigned char *result = new unsigned char[1];
  memset(result, 0x00, 1);
  // int ncell = nrow*ncol;
  // unsigned char *result = new unsigned char[nrow*ncol]; 
  // memset(result, 0x00, ncell);

  /* Read in data */

  int snp_major = start[2];
  int part=0, i=0, j=0;
  size_t ij=0;
  while (1) {
    unsigned char byte;
    if (!part) {
      if (feof(in) || !fread(&byte, 1, 1, in)) {
	printf("Unexpected end of file reached");
	exit(0);
      }
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    // result[ij] = recode[code];
    // returnMat->matrix[i][j] = result[ij];
    result[0] = recode[code];
    returnMat->matrix[i][j] = result[0];
    if(returnMat->matrix[i][j]==3)
      returnMat->matrix[i][j]=1;
    else if(returnMat->matrix[i][j]==1)
      returnMat->matrix[i][j]=3;
    else if(returnMat->matrix[i][j]<0 || returnMat->matrix[i][j]>3){
      printf("Problem in bed file at position=(%d,%d)=%d\n",i,j,returnMat->matrix[i][j]);
      exit(0);
    }
    // printf("(%d,%d)=%d ",i,j,result[ij]);
    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }
  fclose(in);
  delete [] result;
  return returnMat;
}

usiMatrix *bed_to_usiMatrix(const char* file, int nrow,int ncol) {

  //  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  //0,1,2: 3 is missing 
  const unsigned char recode[4] = { '\x02','\x01', '\x03', '\x00'};
  const unsigned char mask = '\x03';


  FILE *in = fopen(file, "r");
  if (!in){
    printf("Couln't open input file: %s", file);
    exit(0);
  }
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3){
    printf("Failed to read first 3 bytes");
    exit(0);
  }
  if (start[0]!='\x6C' || start[1]!='\x1B'){
    printf("Input file does not appear to be a .bed file (%X, %X)", 
	   start[0], start[1]);
    exit(0);
  }
  /* Create output object */
  
  usiMatrix *returnMat = allocUSIntMatrix(nrow,ncol);


  unsigned char *result = new unsigned char[1];
  memset(result, 0x00, 1);
  // int ncell = nrow*ncol;
  // unsigned char *result = new unsigned char[nrow*ncol]; 
  // memset(result, 0x00, ncell);

  /* Read in data */

  int snp_major = start[2];
  int part=0, i=0, j=0;
  size_t ij=0;
  while (1) {
    unsigned char byte;
    if (!part) {
      if (feof(in) || !fread(&byte, 1, 1, in)) {
	printf("Unexpected end of file reached");
	exit(0);
      }
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    // result[ij] = recode[code];
    // returnMat->matrix[i][j] = result[ij];
    result[0] = recode[code];
    returnMat->matrix[i][j] = result[0];
    if(returnMat->matrix[i][j]==3)
      returnMat->matrix[i][j]=1;
    else if(returnMat->matrix[i][j]==1)
      returnMat->matrix[i][j]=3;
    else if(returnMat->matrix[i][j]<0 || returnMat->matrix[i][j]>3){
      printf("Problem in bed file at position=(%d,%d)=%d\n",i,j,returnMat->matrix[i][j]);
      exit(0);
    }
    // printf("(%d,%d)=%d ",i,j,result[ij]);
    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }
  fclose(in);
  delete [] result;
  return returnMat;
}


int numberOfLines(const char *filename){
  const int SIZE=500000;
  char buffer[SIZE];
  int numLines=0;
  std::ifstream pFile (filename,std::ios::in);
  
  if(!pFile){
fileError(filename);
    return 0;
  }
  while(!pFile.eof()){
    pFile.getline(buffer,SIZE);
    numLines++;
  }
  return numLines;
}



/*
we want to cut the first row and the fourth row.
chromosomes and positions.
But we want to check if the cromosomes are 1-22
These are the ones to include, so well return a keeplist,
and the correct pars->chromo and pars->positions.

This is rather slow because we use make a new vector on each line.
This can be optimized in future versions. */
bArray *doBimFile(myPars* pars,const char *filename,const std::string delim){
   ///@param filename A filename to read.@param delim A string of delimiters.
  
  std::vector<int> chromos;//for all lines
  std::vector<double> positions;//for all lines
  

  const int SIZE=500000;//this should only be 6 elements but lets make it big..
  char buffer[SIZE];
 
  std::ifstream pFile (filename,std::ios::in);
 
  
  if(!pFile){
    fileError(filename);
    exit(0);
  }

  
  std::string tmp_string;
  int itemsInRow;
  int numRows =0;
  while(!pFile.eof()){
    pFile.getline(buffer,SIZE);
    tmp_string = std::string(buffer);
    std::vector<std::string> tokens;
    itemsInRow = get_lexemes(tmp_string,tokens,delim);
      
    if (itemsInRow!=6 && itemsInRow!=0){
      printf("plink bim file:%s doesn't have 6 columns in row:%d\n",filename,numRows);
      exit(0);
    }else if(itemsInRow==0)
      break;
    chromos.push_back(atoi((tokens[0]).c_str()));
    positions.push_back(atof((tokens[3]).c_str()));
    numRows++;
  }

  int warn=0;
  //  copy(chromos.begin(), chromos.end(), std::ostream_iterator<int>(std::cout, ", "));  
  //  copy(positions.begin(), positions.end(), std::ostream_iterator<float>(std::cout, ", "));  
  bArray *ret = allocBoolArray(chromos.size());
  int numTrue = 0;
  for(unsigned int i=0;i<chromos.size();i++){
    ret->array[i] = 1;
    numTrue++;
 

  }
  ret->numTrue = numTrue;

  dArray *pos = allocDoubleArray(ret->numTrue);
  iArray *chr = allocIntArray(ret->numTrue);
  int atPos=0;
  for(int i=0;i<ret->x;i++){
    if(ret->array[i]){
      pos->array[atPos] = positions[i]/PLINK_POS_SCALING;
      chr->array[atPos] = chromos[i];
      
      atPos++;
    }

  }
  pars->chr  = chr;
  pars->position= pos;




  return ret;
}



bArray *doUseIndsArray(int nInd, const char *filename, const std::string delim){

  std::ifstream useIndsFile (filename, std::ios::in);

  const int SIZE=2000;
  char buffer[SIZE];
  
  if(!useIndsFile){
    fileError(filename);
    exit(0);
      }

  std::string tmp_string;
  bArray *ret = allocBoolArray(nInd);
  std::vector<int> useIndstmp;
  int numTrue = 0;
  int numRows=0;
  int itemsInRow;
  while(!useIndsFile.eof()){

    useIndsFile.getline(buffer,SIZE);

    tmp_string = std::string(buffer);
    std::vector<std::string> tokens;
    itemsInRow = get_lexemes(tmp_string, tokens, delim);
        if (itemsInRow!=2 && itemsInRow!=0){
     printf("useInds file:%s doesn't have 2 columns in row:%d\n",filename,numRows);
     exit(0);
    }else if(itemsInRow==0)
     break;

    useIndstmp.push_back(atoi(tokens[1].c_str()));
    
    numTrue += useIndstmp[numRows];
    
    numRows++;
    
  }

  for(int i=0;i<nInd;i++)
    ret->array[i] = useIndstmp[i];

  ret-> numTrue = numTrue;

  return ret;
}
