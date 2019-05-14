void getOptions(const char *fileName, functionPars *pars);
dArray *getPos(const char *pName);
iArray *getChr(const char *pName);
iMatrix *getData(const char * pName);
int fexists(std::string str);
std::vector<iArray*> getTestIndividuals(const char *pName);
iMatrix *bed_to_iMatrix(const char* file, int nrow,int ncol) ;
int numberOfLines(const char *filename);
bArray *doBimFile(myPars* pars,const char *filename,const std::string delim,int autosomeMax);
