#include <cmath>
#include "RandomSqrtSingleton.h"


using namespace prompt_gamma_reconstruction;

//Declare statics globally

RandomSqrtSingleton *RandomSqrtSingleton::ptrToSelf = NULL;
RandomSqrtSingleton * RandomSqrtSingleton::Instance(){
    if( 0 == ptrToSelf){
        ptrToSelf = new RandomSqrtSingleton;
    }
    return ptrToSelf;
};

void RandomSqrtSingleton::fillRandomValuesVector(){
    for(size_t i=0; i<RandomSqrtSingleton::NUMBER_OF_POINTS; ++i){
        values.push_back(sqrt(rand_.Rndm()));
    }
    values_vector_end_ = values.end();
    values_vector_begin_ = values.begin();

    itr_current_element_ = values_vector_begin_;
};
