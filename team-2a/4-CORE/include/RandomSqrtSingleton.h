#ifndef RANDOM_SQRT_SINGLETON_H_
#define RANDOM_SQRT_SINGLETON_H_
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>

//Prompt Gamma Includes
#include "ComptonScatter.h"
#include "utilities/Random.h"

namespace prompt_gamma_reconstruction{
    
/*! \brief Returns the sqrt of a random number between 0 and 1.
 * 
 * The random generation and sqrt calculations are performed 
 * during initialization. Users should pay attention to how this
 * class is used because the period is fairly short.
 * 
 * @author Dennis Mackin
 */
class  RandomSqrtSingleton { 

  public:
    inline double getRand(){

      if(++itr_current_element_ == values_vector_end_) itr_current_element_ = values_vector_begin_;

      return (*itr_current_element_);
    };

    static  RandomSqrtSingleton * Instance();

  private:
    pg_tools::Random rand_;
    static const size_t NUMBER_OF_POINTS = 10000103;//arbitrary large prime
    size_t current_value_;
    std::vector<double>::const_iterator itr_current_element_;
    std::vector<double>::const_iterator values_vector_begin_;
    std::vector<double>::const_iterator values_vector_end_;

    std::vector<double> values;

    static RandomSqrtSingleton *ptrToSelf;
    void fillRandomValuesVector();

  protected:
    RandomSqrtSingleton():rand_(12345678),current_value_(0){ fillRandomValuesVector(); };

};//end of class

};//end of namespace
#endif //RANDOM_SQRT_SINGLETON_H_
