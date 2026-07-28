/////////////////////////////////////////////////////////////////
/// Random
/////////////////////////////////////////////////////////////////
/// Wrapper around C++11 random 
///
/// Created 2015-11-03 by Dennis Mackin
/////////////////////////////////////////////////////////////////
#ifndef RANDOM_H_
#define RANDOM_H_

#include <iostream>
#include <random>

using namespace std;
namespace pg_tools {

    template<typename T>
    class Random_T{

      public:
         Random_T(size_t seed=0){
            if(0 == seed){
               random_device r;
               seed = r();
                std::cout <<"WARNING: Random_T called without SEED . . ." << std::endl;
            }
             gen_.seed(seed);
             for(auto i=0; i < 10; ++i) { Rndm();} //warmup the generator
         };

         T Rndm(){
            return generate_canonical<T, std::numeric_limits<T>::digits>(gen_);
         };

        T Gaus(T mean, T sigma){
            std::normal_distribution<T> normal_dist(mean, sigma);
            return normal_dist(gen_);
        };

         T Uniform(T range_min, T range_max){
            std::uniform_real_distribution<> uni(range_min, range_max);
            return uni(gen_);
         };

        void SetSeed(size_t seed){ gen_.seed(seed); };


      private:
         minstd_rand gen_;
    };

    typedef Random_T<float> Random;
};
#endif //RANDOM_H_
