#include "RandomSingleton.h"

using namespace prompt_gamma_reconstruction;

RandomSingleton *RandomSingleton::ptrToSelf = NULL;
RandomSingleton * RandomSingleton::Instance(){
  if( 0 == ptrToSelf){
    ptrToSelf = new RandomSingleton;
  }
  return ptrToSelf;
};

void RandomSingleton::fillRandomValuesVector(){
    values.reserve(NUMBER_OF_RANDS);
    for(size_t i=0; i<RandomSingleton::NUMBER_OF_RANDS; ++i){
        values.push_back(rand_.Rndm());
    }
};
