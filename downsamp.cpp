#include <iostream>
#include <random>
#include <chrono>


 // function to generate random number between 0 and 1, adapted from Serge Rogatch comment in https://stackoverflow.com/questions/9878965/rand-between-0-and-1
int rval(){
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    
    // ready to generate random numbers
      double RandomNumber = unif(rng);
      
      return RandomNumber;
}


void downsampleSites(int **genos, int **F, int nSites, int **newGenos, int **newF,  int nSitesKeep){

  int keepLeft = nSitesKeep;
  int totalLeft = nSites;
  int p = 0;
  int r = 0;
  for(int j = 0; j<nSites; j++){

    p = keepLeft / totalLeft;

    r = rval();


  }

} 
