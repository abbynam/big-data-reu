#ifndef RANDOM_SINGLETON_H_
#define RANDOM_SINGLETON_H_
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <algorithm>

//LIBRARY INCLUDES

#include "utilities/Random.h"

namespace prompt_gamma_reconstruction{
/*! \brief Fast pre-calculated random number generator
 * 
 * Generates and stores a large prime number of randoms 
 * between 0 and 1. Iterates through the list when random numbers are 
 * requested. 
 * 
 * Significantly faster than standard random number classes, <b>but</b>
 * to make it thread safe, threads should call getRandIndex. Threads
 * must manage the index. 
 * 
 * @author Dennis Mackin
 */
class  RandomSingleton {

  public:
    inline double getRand(){
        if(NUMBER_OF_RANDS == current_value_) current_value_ = 0;
        return values[current_value_++];
    };
    inline double getRandIndex(long long  index){
        return values[index % NUMBER_OF_RANDS];
    };
    inline double getLargestIndex(size_t index){
        return static_cast<double>(NUMBER_OF_RANDS);
    };

    inline double getRandGaus(){
        double x1, x2, w;
        do {
            x1 = 2.0 * getRand() - 1.0;
            x2 = 2.0 * getRand() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        return x1 * w;
    };

    static  RandomSingleton * Instance();

    void setSeed(long seed){
        rand_.SetSeed(seed);
        fillRandomValuesVector();
    }

  private:
    pg_tools::Random rand_;
    static const long NUMBER_OF_RANDS = 10000019; //choose large prime
    size_t current_value_;

    std::vector<double> values;

    static RandomSingleton *ptrToSelf;
    void fillRandomValuesVector();

  protected:
    RandomSingleton():rand_(765432),current_value_(0){ fillRandomValuesVector(); };

};//end of class

};//end of namespace
#endif //RANDOM_SINGLETON_H_

